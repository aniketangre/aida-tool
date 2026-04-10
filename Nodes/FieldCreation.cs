using Synera.Core;
using Synera.Core.Graph.Data;
using Synera.Core.Graph.Enums;
using Synera.Core.Implementation.ApplicationService;
using Synera.Core.Implementation.Graph;
using Synera.DataTypes;
using Synera.Kernels;
using Synera.Kernels.DataTypes;
using Synera.Kernels.Fem.Elements;
using Synera.Kernels.Fem.Model;
using Synera.Kernels.Fem.Results;
using Synera.Kernels.Geometry;
using Synera.Kernels.Mesh;
using Synera.Utilities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Runtime.InteropServices;
using Math = System.Math;
using AidaTool.DataTypes;

namespace AidaTool.Nodes
{
    [Guid("9E9D5542-F7FB-4BF0-995E-D218DF2D312E")]
    public sealed class FieldCreation : Node
    {
        public FieldCreation() : base("Create stress field")
        {
            Category = WobisCategories.Aida;
            Subcategory = WobisSubcategories.Field;
            Description = "This node deforms a input mesh based on the fea stress results.";
            GuiPriority = 0;
            Keywords = "deform mesh;stress-guided";

            // Input parameters
            InputParameterManager.AddParameter<IModel>("Input model", "Solved fea model for which we want to get stress results.", ParameterAccess.Item);
            InputParameterManager.AddParameter<SyneraInt>("Iterations", "Iterations determine how far the nodes travel toward their ideal weighted position.", ParameterAccess.Item);
            InputParameterManager.AddParameter<SyneraInt>("Alpha", "Alpha determines how much the stress values influence the node movement relative to the geometric center.", ParameterAccess.Item);

            // Output parameters
            OutputParameterManager.AddParameter<Vector3D>("Max Principle Stress", "Maximum principle stress vectors.", ParameterAccess.List);
            OutputParameterManager.AddParameter<Vector3D>("Min Principle Stress", "Minimum principle stress vectors.", ParameterAccess.List);
            OutputParameterManager.AddParameter<SyneraDouble>("Von-Mises Stress", "Von-Mises stress vectors.", ParameterAccess.List);
            OutputParameterManager.AddParameter<Point3D>("Vertices", "Locations of the the nodes.", ParameterAccess.List);
            OutputParameterManager.AddParameter<Point3D>("Deformed Vertices", "Vertices deformed by vertices.", ParameterAccess.List);
            OutputParameterManager.AddParameter<IMesh>("Mesh", "Mesh.", ParameterAccess.Item);
        }

        public void DeformNodes(List<INode> allNodes, int iterations, double alpha,
                Dictionary<int, double> nodeAverageStress,
                Dictionary<int, List<int>> nodeToNodesMap,
                HashSet<int> boundaryNodeIndicesSet)
        {
            // 1. Find Min/Max Stress for normalization
            double minStress = double.MaxValue;
            double maxStress = double.MinValue;
            foreach (var node in allNodes)
            {
                minStress = Math.Min(minStress, nodeAverageStress[node.Index]);
                maxStress = Math.Max(maxStress, nodeAverageStress[node.Index]);
            }

            // Create node dictionary
            Dictionary<int, INode> nodesDict = new Dictionary<int, INode>(allNodes.Count);
            foreach (INode node in allNodes)
            {
                nodesDict[node.Index] = node;
            }

            // 2. Iterative Relaxation
            for (int iter = 0; iter < iterations; iter++)
            {
                for (int i = 0; i < allNodes.Count; i++)
                {
                    // Keep boundary nodes fixed
                    if (boundaryNodeIndicesSet.Contains(allNodes[i].Index))
                    {
                        allNodes[i] = allNodes[i];
                        continue;
                    }

                    Point3D weightedSum = new Point3D(0, 0, 0);
                    double totalWeight = 0;
                    foreach (var connectedNodeIdx in nodeToNodesMap[allNodes[i].Index])
                    {
                        // Calculate weight based on neighbor's stress
                        double normalizedStress = (nodeAverageStress[connectedNodeIdx] - minStress) / (maxStress - minStress + 1e-6f);
                        double weight = 1.0f + (alpha * normalizedStress);

                        weightedSum += weight * nodesDict[connectedNodeIdx].Location;
                        totalWeight += weight;
                    }

                    allNodes[i].Location = weightedSum / totalWeight;
                }
            }
        }

        // Function to calculate von mises stress
        public double CalculateVonMises(double[] stressTensor)
        {
            if (stressTensor == null)
            {
                throw new ArgumentNullException(nameof(stressTensor));
            }
            if (stressTensor.Length < 6)
            {
                throw new ArgumentException("Tensor array must contain at least 6 components: [Sxx, Syy, Szz, Sxy, Syz, Sxz].", nameof(stressTensor));
            }

            double sxx = stressTensor[0];
            double syy = stressTensor[1];
            double szz = stressTensor[2];
            double sxy = stressTensor[3];
            double syz = stressTensor[4];
            double sxz = stressTensor[5];

            double normalPart = (sxx - syy) * (sxx - syy) + (syy - szz) * (syy - szz) + (szz - sxx) * (szz - sxx);
            double shearPart = 6.0 * ((sxy * sxy) + (syz * syz) + (sxz * sxz));

            return Math.Sqrt(0.5 * (normalPart + shearPart));
        }

        protected override void SolveInstance(IDataAccess dataAccess)
        {
            bool isDataSuccess = dataAccess.GetData(0, out IModel model);
            isDataSuccess &= dataAccess.GetData(1, out SyneraInt iterations);
            isDataSuccess &= dataAccess.GetData(2, out SyneraInt alpha);

            if (!isDataSuccess)
                return;

            // Validation of inputs
            bool flag = true;
            // Check input fem model
            if (!model.IsValid)
            {
                AddError(0, "Invalid model.");
                flag &= false;
            }
            else if (!(model.HasShellElementsOnly()))
            {
                AddError(0, "Model with only shell elements is expected.");
                flag &= false;
            }
            // validation for iterations
            if (iterations <= 0)
            {
                AddError(1, "Iterations must be a positive integer.");
                flag &= false;
            }
            // validation for alpha
            if (alpha < 0)
            {
                AddError(2, "Alpha must be a non-negative integer.");
                flag &= false;
            }
            if (double.IsPositiveInfinity(alpha))
            {
                AddError(2, "Node failed because given given alpha value is infinite. Try to input alpha such that it is valid number");
                flag &= false;
            }

            if (!flag)
                return;

            // Function to calculate von mises stress
            //static double CalculateVonMises(double[] stressTensor)
            //{
            //    if (stressTensor == null)
            //    {
            //        throw new ArgumentNullException(nameof(stressTensor));
            //    }
            //    if (stressTensor.Length < 6)
            //    {
            //        throw new ArgumentException("Tensor array must contain at least 6 components: [Sxx, Syy, Szz, Sxy, Syz, Sxz].", nameof(stressTensor));
            //    }

            //    double sxx = stressTensor[0];
            //    double syy = stressTensor[1];
            //    double szz = stressTensor[2];
            //    double sxy = stressTensor[3];
            //    double syz = stressTensor[4];
            //    double sxz = stressTensor[5];

            //    double normalPart = (sxx - syy) * (sxx - syy) + (syy - szz) * (syy - szz) + (szz - sxx) * (szz - sxx);
            //    double shearPart = 6.0 * ((sxy * sxy) + (syz * syz) + (sxz * sxz));

            //    return Math.Sqrt(0.5 * (normalPart + shearPart));
            //}

            // Function to Refine mesh
            //static void DeformNodes(Dictionary<int, INode> allNodes, int iterations, double alpha,
            //    Dictionary<int, double> nodeAverageStress,
            //    Dictionary<int, List<int>> nodeToNodesMap,
            //    HashSet<int> boundaryNodeIndicesSet)
            //{
            //    // 1. Find Min/Max Stress for normalization
            //    double minStress = double.MaxValue;
            //    double maxStress = double.MinValue;

            //    foreach (var (idx, node) in allNodes)
            //    {
            //        minStress = Math.Min(minStress, nodeAverageStress[idx]);
            //        maxStress = Math.Max(maxStress, nodeAverageStress[idx]);
            //    }

            //    // 2. Iterative Relaxation
            //    for (int iter = 0; iter < iterations; iter++)
            //    {
            //        Dictionary<int, Point3D> newPositions = new Dictionary<int, Point3D>(allNodes.Count);

            //        foreach (var (idx, node) in allNodes)
            //        {
            //            // Keep boundary nodes fixed
            //            if (boundaryNodeIndicesSet.Contains(node.Index))
            //            {
            //                newPositions[idx] = node.Location;
            //                continue;
            //            }

            //            Point3D weightedSum = new Point3D(0, 0, 0);
            //            double totalWeight = 0;
            //            foreach (var connectedNodeIdx in nodeToNodesMap[node.Index])
            //            {
            //                // Calculate weight based on neighbor's stress
            //                double normalizedStress = (nodeAverageStress[connectedNodeIdx] - minStress) / (maxStress - minStress + 1e-6f);
            //                double weight = 1.0f + (alpha * normalizedStress);

            //                weightedSum += weight * allNodes[connectedNodeIdx].Location;
            //                totalWeight += weight;
            //            }
            //            newPositions[idx] = weightedSum / totalWeight;
            //        }

            //        // 3. Update positions after calculating all new points (Synchronous update)
            //        foreach (var (idx, node) in allNodes)
            //        {
            //            allNodes[idx].Location = newPositions[idx];
            //        }
            //    }
            //}

            // Initialize progress tracking
            Progress operationalProgressMeter = this.CreateProgress();
            StressStreamlineBuilder StressStreamlineBuilder = new StressStreamlineBuilder(model);

            // Calculate principal stress directions (40% of progress)
            Progress subTaskOneProgress = operationalProgressMeter.CreateSubtask(0.5);
            Dictionary<int, (Vector3D majorVector, Vector3D minorVector)> stressVectors = StressStreamlineBuilder.GetPrincipalStressValues(subTaskOneProgress);

            // Extract the major and minor principal stress vectors into separate lists
            List<Vector3D> majorPrincipalStressVectors = stressVectors.Values.Select(v => v.majorVector).ToList();
            List<Vector3D> minorPrincipalStressVectors = stressVectors.Values.Select(v => v.minorVector).ToList();

            // Get all shell elements and nodes from the model
            IEnumerable<IShellElement> shellElements = model.Elements.OfType<IShellElement>();
            List<INode> nodes = model.Nodes.ToList();

            // Get a dictionary of node id and ids of the elements connected to it
            Dictionary<int, List<int>> nodeToElementsMap = new Dictionary<int, List<int>>();
            foreach (IShellElement shellElement in shellElements)
            {
                List<int> cornerNodeIndices = shellElement.NodeIndices.Except(shellElement.GetMidNodeIndices()).ToList();

                // Write me a code to fill the nodeToElementsMap dictionary with the node indices and the corresponding element indices
                foreach (int nodeIndex in cornerNodeIndices)
                {
                    if (!nodeToElementsMap.ContainsKey(nodeIndex))
                    {
                        nodeToElementsMap[nodeIndex] = new List<int>();
                    }
                    nodeToElementsMap[nodeIndex].Add(shellElement.Index);
                }
            }

            // Get a dictionary of node id and the ids of it's neighboring nodes.
            Dictionary<int, List<int>> nodeToNeighboringNodesMap = new Dictionary<int, List<int>>();
            foreach (IndexPair edge in model.Edges)
            {
                int a = edge.A;
                int b = edge.B;

                // skip degenerate edges
                if (a == b)
                    continue;

                if (!nodeToNeighboringNodesMap.ContainsKey(a))
                    nodeToNeighboringNodesMap[a] = new List<int>();
                if (!nodeToNeighboringNodesMap.ContainsKey(b))
                    nodeToNeighboringNodesMap[b] = new List<int>();

                // add neighbor (prevent duplicates)
                if (!nodeToNeighboringNodesMap[a].Contains(b))
                    nodeToNeighboringNodesMap[a].Add(b);
                if (!nodeToNeighboringNodesMap[b].Contains(a))
                    nodeToNeighboringNodesMap[b].Add(a);
            }

            // Get a dictionary of node id and the ids of it's neighboring nodes.
            //Dictionary<int, HashSet<int>> nodeToNeighboringNodesMap = new Dictionary<int, HashSet<int>>(nodes.Count);
            //foreach (IShellElement shellElement in shellElements)
            //{
            //    List<int> cornerNodeIndices = shellElement.NodeIndices.Except(shellElement.GetMidNodeIndices()).ToList();

            //    for (int i = 0; i < cornerNodeIndices.Count; ++i)
            //    {
            //        if (!nodeToNeighboringNodesMap.TryGetValue(cornerNodeIndices[i], out var connectedNodes))
            //        {
            //            connectedNodes = new HashSet<int>();
            //            nodeToNeighboringNodesMap[cornerNodeIndices[i]] = connectedNodes;
            //        }

            //        // add all other corner nodes of the element as neighbors (skip self)
            //        for (int j = 0; j < cornerNodeIndices.Count; ++j)
            //        {
            //            if (i == j) continue;
            //            connectedNodes.Add(cornerNodeIndices[j]);
            //        }
            //    }
            //}

            // Get a element stress results for the first time step (assuming results are available and properly indexed)
            IResultsAtTimeStep elementStressResults = model.Results
                .First(result => result.Name.Equals(ResultType.ElementStressTensor.Name))
                .Values.First() // Get the first result set (e.g., for a single load case)
                .Values.First(); // Get the results for the first time step

            // Get a dictionary of average stress values for each node based on the elements connected to it
            Dictionary<int, double> nodeAverageStress = new Dictionary<int, double>(nodes.Count);
            foreach (INode node in nodes)
            {
                if (nodeToElementsMap.ContainsKey(node.Index))
                {
                    List<int> connectedElementIndices = nodeToElementsMap[node.Index];
                    double sumStress = 0;
                    foreach (int idx in connectedElementIndices)
                    {
                        sumStress += CalculateVonMises(elementStressResults[idx]);
                    }
                    nodeAverageStress[node.Index] = sumStress / connectedElementIndices.Count;
                }
                else
                {
                    nodeAverageStress[node.Index] = 0; // If no connected elements, set average stress to 0
                }
            }

            // Get a hash set of boundary nodes
            IEnumerable<int> nakedEdgeIndices = StressStreamlineBuilder.Mesh.GetNakedEdgeIndices();
            HashSet<int> boundaryNodeIndicesSet = new HashSet<int>();
            foreach (int idx in nakedEdgeIndices)
            {
                IndexPair edge = model.GetNodeIndicesForEdge(idx);
                boundaryNodeIndicesSet.AddRange(edge.ToList());
            }
            List<Point3D> boundaryNodes = (from node in nodes
                                          where boundaryNodeIndicesSet.Contains(node.Index)
                                          select node.Location).ToList();

            // Calculate element centers for all elements.
            List<Point3D> elementCenters = shellElements.Select(elem => elem.GetCenter()).ToList();

            // Get node locations for the output (optional, can be used for visualization or further processing)
            List<Point3D> nodeLocations = nodes.Select(node => node.Location).ToList();

            // Refined mesh
            List<INode> originalNodes = nodes;
            DeformNodes(originalNodes, iterations, alpha, nodeAverageStress, nodeToNeighboringNodesMap, boundaryNodeIndicesSet);
            List<Point3D> updatedNodes = originalNodes.Select(x => x.Location).ToList();
            List<Point3D> deformedVertices = originalNodes
                                        .Select(node => node.Location) // Extract the Location
                                        .ToList();

            // Create mesh
            List<MeshFace> meshFaces = new List<MeshFace>(shellElements.Count());
            Dictionary<int, IShellElement> shellElementsDict = new Dictionary<int, IShellElement>(shellElements.Count());
            foreach (IShellElement shellElement in shellElements)
            {
                shellElementsDict[shellElement.Index] = shellElement;
            }

            foreach (var shellElement in shellElements)
            {
                meshFaces.Add(new MeshFace(shellElement.NodeIndices.Select(x => x-1).ToArray()));
            }

            // Create a mesh 
            IMeshKernel meshKernel = Application.Current.KernelManager.Get<IMeshKernel>();
            IMesh mesh = meshKernel.CreateFromVerticesAndFaces(deformedVertices, meshFaces);

            // Set the output data
            dataAccess.SetListData(0, majorPrincipalStressVectors);
            dataAccess.SetListData(1, minorPrincipalStressVectors);
            dataAccess.SetListData(2, nodeAverageStress.Values.ToList());
            dataAccess.SetListData(3, nodeLocations);
            dataAccess.SetListData(4, deformedVertices);
            dataAccess.SetData(5, mesh);
        }
    }
}

