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
using Synera.Kernels.Mesh;
using Synera.Utilities;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using Math = System.Math;
using AidaTool.DataTypes;

namespace AidaTool.Nodes
{
    // Each Synera node needs a unique GUID so the software can identify it persistently.
    // If you change this GUID, existing graphs that use this node will lose the connection.
    [Guid("3D8F2A1B-5C6E-4D7F-8E9A-0B1C2D3E4F50")]
    public sealed class QuadMeshFD : Node
    {
        // The constructor sets up everything the user sees in the Synera UI:
        // the node name, category, description, and all input/output ports.
        public QuadMeshFD() : base("Quad Mesh FD")
        {
            Category = WobisCategories.Aida;
            Subcategory = WobisSubcategories.Field;
            Description = "Deforms a mesh guided by both von Mises stress magnitude and principal stress " +
                          "directions. Inspired by graphic statics, the relaxation is anisotropic: nodes " +
                          "are attracted more strongly toward neighbors that lie along the principal stress " +
                          "trajectories, causing mesh elements to elongate in those directions.";
            GuiPriority = 1;
            Keywords = "deform mesh;stress-guided;principal stress;graphic statics;anisotropic";

            // ── Input ports ───────────────────────────────────────────────────────────────
            // Port index 0: the solved FEA model containing geometry and stress results.
            InputParameterManager.AddParameter<IModel>("Input model",
                "Solved FEA model whose stress results will guide the mesh deformation.",
                ParameterAccess.Item);

            // Port index 1: how many times the relaxation loop runs.
            // More iterations → smoother result, but slower to compute.
            InputParameterManager.AddParameter<SyneraInt>("Iterations",
                "Number of relaxation iterations. Higher values produce smoother, more stress-aligned meshes.",
                ParameterAccess.Item);

            // Port index 2: controls how much von Mises stress pulls nodes toward high-stress areas.
            // Alpha = 0 means all neighbors are weighted equally (standard Laplacian smoothing).
            InputParameterManager.AddParameter<SyneraInt>("Alpha",
                "Controls how strongly von Mises stress influences node movement relative to the geometric center. " +
                "Alpha = 0 gives uniform Laplacian smoothing; higher values pull nodes toward high-stress regions.",
                ParameterAccess.Item);

            // Port index 3: controls how much the principal stress directions steer node movement.
            // Beta = 0 means no directional effect — identical output to FieldCreation.
            // Higher Beta rotates the displacement direction toward stress trajectories
            // without changing how far the node moves (no overshoot risk).
            InputParameterManager.AddParameter<SyneraDouble>("Beta",
                "Controls how strongly principal stress directions influence node movement. " +
                "Beta = 0 disables directional anisotropy; higher values rotate the displacement " +
                "progressively toward the principal stress axes without changing the step size, " +
                "so mesh inversion is not possible.",
                ParameterAccess.Item);

            // ── Output ports ──────────────────────────────────────────────────────────────
            // These expose the intermediate stress data and the final deformed mesh.

            // Port index 0: one vector per element pointing along the major principal stress axis.
            OutputParameterManager.AddParameter<Vector3D>("Max Principle Stress",
                "Maximum principal stress vectors per element.",
                ParameterAccess.List);

            // Port index 1: one vector per element pointing along the minor principal stress axis.
            OutputParameterManager.AddParameter<Vector3D>("Min Principle Stress",
                "Minimum principal stress vectors per element.",
                ParameterAccess.List);

            // Port index 2: von Mises stress scalar averaged to each node from its surrounding elements.
            OutputParameterManager.AddParameter<SyneraDouble>("Von-Mises Stress",
                "Von Mises stress value averaged per node.",
                ParameterAccess.List);

            // Port index 3: node positions before any deformation — useful for comparison.
            OutputParameterManager.AddParameter<Point3D>("Vertices",
                "Original node locations before deformation.",
                ParameterAccess.List);

            // Port index 4: node positions after the stress-guided relaxation.
            OutputParameterManager.AddParameter<Point3D>("Deformed Vertices",
                "Node locations after stress-guided anisotropic deformation.",
                ParameterAccess.List);

            // Port index 5: the final deformed mesh built from the updated node positions.
            OutputParameterManager.AddParameter<IMesh>("Mesh",
                "Resulting deformed mesh aligned with the principal stress field.",
                ParameterAccess.Item);
        }

        /// <summary>
        /// Converts a 6-component stress tensor into a single von Mises stress scalar.
        ///
        /// The von Mises criterion combines all stress components into one value that
        /// represents the overall "intensity" of the stress state at a point. It is widely
        /// used in structural engineering to predict yielding.
        ///
        /// The formula separates the tensor into:
        ///   - Normal stress differences (how much the material is being stretched unevenly)
        ///   - Shear stresses (how much the material is being twisted)
        ///
        /// Input order: [Sxx, Syy, Szz, Sxy, Syz, Sxz]
        /// </summary>
        private static double CalculateVonMisesStress(double[] stressTensor)
        {
            // Extract individual stress components from the tensor array.
            double sxx = stressTensor[0], syy = stressTensor[1], szz = stressTensor[2];
            double sxy = stressTensor[3], syz = stressTensor[4], sxz = stressTensor[5];

            // Normal stress contribution: differences between axial stress components.
            double normalPart = (sxx - syy) * (sxx - syy)
                              + (syy - szz) * (syy - szz)
                              + (szz - sxx) * (szz - sxx);

            // Shear stress contribution: all three shear components combined.
            double shearPart = 6.0 * (sxy * sxy + syz * syz + sxz * sxz);

            // Final von Mises stress — always a positive scalar.
            return Math.Sqrt(0.5 * (normalPart + shearPart));
        }

        /// <summary>
        /// Moves interior mesh nodes iteratively so the mesh aligns with the stress field.
        /// Uses a Gauss-Seidel update: each node's new position is applied immediately,
        /// so the next node in the loop already sees the updated neighbor positions.
        /// This converges faster than Jacobi and is stable for quad meshes.
        ///
        /// Each node moves in two stages:
        ///   Stage 1 — Von Mises weighted centroid:
        ///     Compute a weighted average position of the node's neighbors.
        ///     Neighbors in high-stress regions get a higher weight (controlled by Alpha),
        ///     so the node is pulled toward stress concentrations.
        ///
        ///   Stage 2 — Anisotropic direction steering:
        ///     Decompose the displacement (node → centroid) into components along the
        ///     major and minor principal stress directions and a residual component.
        ///     Amplify the stress-aligned components by Beta, then rescale the entire
        ///     vector back to the original length. This steers the movement direction
        ///     toward the stress trajectories without increasing the step size.
        /// </summary>
        private static void DeformNodes(
            List<INode> allNodes,           // All nodes in the mesh
            int iterations,                 // How many relaxation passes to run
            double alpha,                   // Strength of von Mises stress influence
            double beta,                    // Strength of principal stress direction influence
            Dictionary<int, double> nodeVonMisesStress,             // Von Mises stress per node index
            Dictionary<int, List<int>> nodeToNeighborsMap,          // Neighbor node indices per node
            HashSet<int> boundaryNodeIndices,                        // Nodes that must not move
            Dictionary<int, (Vector3D majorDirection, Vector3D minorDirection)> nodePrincipalDirections) // Stress directions per node
        {
            // Normalize von Mises stress to [0, 1] so Alpha behaves consistently
            // regardless of the absolute stress magnitudes in the model.
            double minStress = nodeVonMisesStress.Values.Min();
            double maxStress = nodeVonMisesStress.Values.Max();
            double stressRange = maxStress - minStress + 1e-6; // small epsilon prevents division by zero

            // Build a dictionary for fast node lookup by index during the inner loop.
            // This avoids repeatedly scanning the list to find a node by its index.
            Dictionary<int, INode> nodeDict = allNodes.ToDictionary(n => n.Index);

            for (int iteration = 0; iteration < iterations; iteration++)
            {
                foreach (INode node in allNodes)
                {
                    // Boundary nodes (on free edges) stay fixed to preserve the mesh outline.
                    if (boundaryNodeIndices.Contains(node.Index))
                        continue;

                    List<int> neighbors = nodeToNeighborsMap[node.Index];
                    if (neighbors.Count == 0)
                        continue;

                    // Retrieve the averaged principal stress directions for this node.
                    // These were computed by averaging over all elements connected to the node.
                    nodePrincipalDirections.TryGetValue(node.Index, out var principalDirs);
                    Vector3D majorDir = principalDirs.majorDirection;
                    Vector3D minorDir = principalDirs.minorDirection;

                    double majorMagnitude = majorDir.Length;
                    double minorMagnitude = minorDir.Length;

                    // Normalize principal magnitudes per node using the local maximum.
                    // This ensures the dominant stress direction always contributes fully (score = 1.0)
                    // and the other direction contributes proportionally.
                    // Using a global maximum would cause stress peaks to suppress all other nodes.
                    double localMax = Math.Max(majorMagnitude, minorMagnitude) + 1e-6;
                    double normalizedMajorMagnitude = majorMagnitude / localMax;
                    double normalizedMinorMagnitude = minorMagnitude / localMax;

                    // Compute unit vectors along each principal stress direction.
                    // If a node has negligible stress (e.g., in a dead zone), fall back to
                    // world-axis directions so the math doesn't break down.
                    Vector3D majorUnit = majorMagnitude > 1e-10 ? majorDir.Normalized() : new Vector3D(1, 0, 0);
                    Vector3D minorUnit = minorMagnitude > 1e-10 ? minorDir.Normalized() : new Vector3D(0, 1, 0);

                    // ── Stage 1: Compute von Mises weighted Laplacian centroid ────────────────
                    // The centroid is a weighted average of neighbor positions.
                    // Neighbors with higher von Mises stress get a larger weight,
                    // pulling the current node toward high-stress regions.
                    Point3D weightedPositionSum = new Point3D(0, 0, 0);
                    double totalWeight = 0.0;

                    foreach (int neighborIdx in neighbors)
                    {
                        // Gauss-Seidel: read the live position from nodeDict.
                        // If this neighbor was already updated earlier in this iteration,
                        // we use its new position — this is what makes Gauss-Seidel faster.
                        INode neighbor = nodeDict[neighborIdx];
                        double normalizedStress = (nodeVonMisesStress[neighborIdx] - minStress) / stressRange;

                        // Weight formula: baseline of 1.0 plus Alpha-scaled stress influence.
                        // A neighbor with maximum stress gets weight = 1 + Alpha.
                        // A neighbor with minimum stress gets weight = 1.0.
                        double vonMisesWeight = 1.0 + alpha * normalizedStress;

                        weightedPositionSum += vonMisesWeight * neighbor.Location;
                        totalWeight += vonMisesWeight;
                    }

                    if (totalWeight <= 1e-10)
                        continue;

                    // The centroid is where the node would move under pure von Mises weighted smoothing.
                    Point3D centroid = weightedPositionSum / totalWeight;

                    // ── Stage 2: Steer displacement direction toward principal stress axes ─────
                    // Instead of moving straight to the centroid, we tilt the movement direction
                    // to align more with the principal stress trajectories.

                    // The displacement vector from current position to centroid.
                    Vector3D delta = Point3D.Subtract(centroid, node.Location);

                    // Project delta onto each principal stress direction to get scalar components.
                    // dMajor > 0 means the centroid is "ahead" along the major stress axis.
                    double dMajor = delta * majorUnit;
                    double dMinor = delta * minorUnit;

                    // Whatever part of delta is not along major or minor stress axes.
                    // For shell elements this is mostly the through-thickness component.
                    Vector3D deltaRest = delta - dMajor * majorUnit - dMinor * minorUnit;

                    // Amplify the stress-aligned components by Beta.
                    // A higher Beta means the node steers more aggressively toward stress directions.
                    // The factor (1 + Beta * normalizedMagnitude) scales each component:
                    //   - Beta = 0: no amplification, direction unchanged
                    //   - Beta > 0: stress-aligned components grow relative to deltaRest
                    Vector3D anisotropicDirection = deltaRest
                        + (1.0 + beta * normalizedMajorMagnitude) * dMajor * majorUnit
                        + (1.0 + beta * normalizedMinorMagnitude) * dMinor * minorUnit;

                    // Rescale the anisotropic direction back to the original delta length.
                    // This is critical: Beta changes DIRECTION only, not step SIZE.
                    // Without this rescaling, high Beta would move nodes further than their centroid,
                    // causing overshoot and potential mesh inversion.
                    double deltaLength = delta.Length;
                    double anisotropicLength = anisotropicDirection.Length;

                    if (anisotropicLength > 1e-10 && deltaLength > 1e-10)
                        node.Location = node.Location + (deltaLength / anisotropicLength) * anisotropicDirection;
                    else
                        node.Location = node.Location + delta; // fallback: move straight to centroid
                }
            }
        }

        /// <summary>
        /// Called by Synera each time the node needs to recompute its outputs.
        /// Reads inputs, validates them, runs the stress-guided relaxation, and sets outputs.
        /// </summary>
        protected override void SolveInstance(IDataAccess dataAccess)
        {
            // ── Read inputs from the connected ports ──────────────────────────────────────
            bool isDataSuccess = dataAccess.GetData(0, out IModel model);
            isDataSuccess &= dataAccess.GetData(1, out SyneraInt iterations);
            isDataSuccess &= dataAccess.GetData(2, out SyneraInt alpha);
            isDataSuccess &= dataAccess.GetData(3, out double beta);

            // If any required input is missing or disconnected, stop here.
            if (!isDataSuccess)
                return;

            // ── Validate inputs ───────────────────────────────────────────────────────────
            // We collect all errors before returning so the user sees all problems at once.
            bool isValid = true;

            if (!model.IsValid)
            {
                AddError(0, "Invalid model.");
                isValid = false;
            }
            else if (!model.HasShellElementsOnly())
            {
                // This node is designed for shell (surface) elements only.
                // Solid elements have a different stress distribution and topology.
                AddError(0, "A model containing shell elements only is expected.");
                isValid = false;
            }

            if (iterations <= 0)
            {
                AddError(1, "Iterations must be a positive integer.");
                isValid = false;
            }

            if (alpha < 0)
            {
                AddError(2, "Alpha must be a non-negative integer.");
                isValid = false;
            }

            if (double.IsPositiveInfinity(alpha))
            {
                AddError(2, "Alpha is infinite. Please provide a finite value.");
                isValid = false;
            }

            if (beta < 0)
            {
                AddError(3, "Beta must be a non-negative value.");
                isValid = false;
            }

            if (double.IsPositiveInfinity(beta))
            {
                AddError(3, "Beta is infinite. Please provide a finite value.");
                isValid = false;
            }

            if (!isValid)
                return;

            // ── Step 1: Compute principal stress vectors per element ──────────────────────
            // StressStreamlineBuilder handles the eigenvalue decomposition of the stress tensor
            // to extract the two in-plane principal stress directions for each shell element.
            Progress progress = this.CreateProgress();
            StressStreamlineBuilder streamlineBuilder = new StressStreamlineBuilder(model);

            // This computation takes roughly half the total processing time.
            Progress stressProgress = progress.CreateSubtask(0.5);
            Dictionary<int, (Vector3D majorVector, Vector3D minorVector)> elementPrincipalStress =
                streamlineBuilder.GetPrincipalStressValues(stressProgress);

            // Separate into flat lists for the output ports.
            List<Vector3D> majorStressVectors = elementPrincipalStress.Values.Select(v => v.majorVector).ToList();
            List<Vector3D> minorStressVectors  = elementPrincipalStress.Values.Select(v => v.minorVector).ToList();

            // ── Step 2: Build mesh topology lookup maps ───────────────────────────────────
            IEnumerable<IShellElement> shellElements = model.Elements.OfType<IShellElement>();
            List<INode> nodes = model.Nodes.ToList();

            // Map: node index → list of element indices that share this node.
            // Needed to average element-level stress data down to each node.
            // Mid-side nodes of higher-order elements (e.g. 8-node quads) are excluded
            // because we only want the corner nodes that define the element shape.
            Dictionary<int, List<int>> nodeToElementsMap = new Dictionary<int, List<int>>();
            foreach (IShellElement element in shellElements)
            {
                List<int> cornerNodeIndices = element.NodeIndices.Except(element.GetMidNodeIndices()).ToList();
                foreach (int nodeIndex in cornerNodeIndices)
                {
                    if (!nodeToElementsMap.ContainsKey(nodeIndex))
                        nodeToElementsMap[nodeIndex] = new List<int>();
                    nodeToElementsMap[nodeIndex].Add(element.Index);
                }
            }

            // Map: node index → list of directly connected neighbor node indices.
            // Two nodes are neighbors if they share a mesh edge.
            // This defines the "neighborhood" used in the Laplacian relaxation.
            Dictionary<int, List<int>> nodeToNeighborsMap = new Dictionary<int, List<int>>();
            foreach (IndexPair edge in model.Edges)
            {
                int a = edge.A;
                int b = edge.B;

                // Skip degenerate edges where both endpoints are the same node.
                if (a == b)
                    continue;

                if (!nodeToNeighborsMap.ContainsKey(a)) nodeToNeighborsMap[a] = new List<int>();
                if (!nodeToNeighborsMap.ContainsKey(b)) nodeToNeighborsMap[b] = new List<int>();

                // Add each node as a neighbor of the other (undirected edge).
                // Check for duplicates to avoid counting shared edges twice.
                if (!nodeToNeighborsMap[a].Contains(b)) nodeToNeighborsMap[a].Add(b);
                if (!nodeToNeighborsMap[b].Contains(a)) nodeToNeighborsMap[b].Add(a);
            }

            // ── Step 3: Compute per-node von Mises stress ─────────────────────────────────
            // FEA stress results are stored per element (at integration points).
            // We average the von Mises stress of all elements connected to each node
            // to get a smooth scalar field defined at the nodes.
            IResultsAtTimeStep elementStressResults = model.Results
                .First(r => r.Name.Equals(ResultType.ElementStressTensor.Name))
                .Values.First()  // first load case
                .Values.First(); // first time step

            Dictionary<int, double> nodeVonMisesStress = new Dictionary<int, double>(nodes.Count);
            foreach (INode node in nodes)
            {
                if (nodeToElementsMap.TryGetValue(node.Index, out List<int> connectedElements))
                {
                    // Sum up von Mises stress from all connected elements, then average.
                    double stressSum = connectedElements
                        .Sum(elemIdx => CalculateVonMisesStress(elementStressResults[elemIdx]));
                    nodeVonMisesStress[node.Index] = stressSum / connectedElements.Count;
                }
                else
                {
                    // Isolated nodes with no connected elements carry zero stress.
                    nodeVonMisesStress[node.Index] = 0.0;
                }
            }

            // ── Step 4: Average principal stress directions per node ───────────────────────
            // Similarly, the element-level principal stress vectors are averaged over all
            // elements connected to each node, producing a smooth directional field at the nodes.
            Dictionary<int, (Vector3D majorDirection, Vector3D minorDirection)> nodePrincipalDirections =
                new Dictionary<int, (Vector3D, Vector3D)>(nodes.Count);

            foreach (INode node in nodes)
            {
                if (!nodeToElementsMap.TryGetValue(node.Index, out List<int> connectedElements))
                    continue;

                Vector3D majorSum = new Vector3D(0, 0, 0);
                Vector3D minorSum = new Vector3D(0, 0, 0);
                int contributingElements = 0;

                foreach (int elemIdx in connectedElements)
                {
                    if (elementPrincipalStress.TryGetValue(elemIdx, out var stressVec))
                    {
                        majorSum = majorSum + stressVec.majorVector;
                        minorSum = minorSum + stressVec.minorVector;
                        contributingElements++;
                    }
                }

                if (contributingElements > 0)
                {
                    // Divide by the count to get the average direction vector.
                    double inv = 1.0 / contributingElements;
                    nodePrincipalDirections[node.Index] = (inv * majorSum, inv * minorSum);
                }
            }

            // ── Step 5: Identify boundary nodes that must not move ────────────────────────
            // Nodes on "naked" (free/open) edges sit on the boundary of the mesh.
            // Keeping them fixed preserves the overall shape and footprint of the structure.
            IEnumerable<int> nakedEdgeIndices = streamlineBuilder.Mesh.GetNakedEdgeIndices();
            HashSet<int> boundaryNodeIndices = new HashSet<int>();
            foreach (int edgeIdx in nakedEdgeIndices)
            {
                IndexPair edge = model.GetNodeIndicesForEdge(edgeIdx);
                boundaryNodeIndices.AddRange(edge.ToList());
            }

            // Record original node positions before deformation for the output port.
            List<Point3D> originalVertices = nodes.Select(n => n.Location).ToList();

            // ── Step 6: Run the stress-guided relaxation ──────────────────────────────────
            // This modifies node positions in-place (Gauss-Seidel update).
            DeformNodes(
                nodes,
                iterations,
                alpha,
                beta,
                nodeVonMisesStress,
                nodeToNeighborsMap,
                boundaryNodeIndices,
                nodePrincipalDirections);

            // Collect updated node positions after deformation.
            List<Point3D> deformedVertices = nodes.Select(n => n.Location).ToList();

            // ── Step 7: Rebuild the mesh from deformed vertex positions ───────────────────
            // The mesh topology (which nodes form which faces) stays the same.
            // Only the vertex positions change.
            List<MeshFace> meshFaces = shellElements
                .Select(elem => new MeshFace(elem.NodeIndices.Select(i => i - 1).ToArray()))
                .ToList();

            IMeshKernel meshKernel = Application.Current.KernelManager.Get<IMeshKernel>();
            IMesh deformedMesh = meshKernel.CreateFromVerticesAndFaces(deformedVertices, meshFaces);

            // ── Step 8: Set all output ports ──────────────────────────────────────────────
            dataAccess.SetListData(0, majorStressVectors);
            dataAccess.SetListData(1, minorStressVectors);
            dataAccess.SetListData(2, nodeVonMisesStress.Values.ToList());
            dataAccess.SetListData(3, originalVertices);
            dataAccess.SetListData(4, deformedVertices);
            dataAccess.SetData(5, deformedMesh);
        }
    }
}
