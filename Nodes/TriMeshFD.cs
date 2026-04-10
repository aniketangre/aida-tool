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
    // Unique identifier for this node. Never change this once the node is in use —
    // Synera uses it to reconnect saved graphs to the correct node class.
    [Guid("4A9B3C2D-6E7F-5A8B-9C0D-1E2F3A4B5C6D")]
    public sealed class TriMeshFD : Node
    {
        // The constructor defines everything the user sees: node name, description,
        // input ports, and output ports.
        public TriMeshFD() : base("Tri Mesh FD")
        {
            Category = WobisCategories.Aida;
            Subcategory = WobisSubcategories.Field;
            Description = "Deforms a triangular shell mesh guided by both von Mises stress magnitude and " +
                          "principal stress directions. Uses a synchronous Jacobi update with per-iteration " +
                          "damping to prevent the element distortion and chaotic node movement that occur " +
                          "when anisotropic relaxation is applied to triangular topology.";
            GuiPriority = 1;
            Keywords = "deform mesh;stress-guided;principal stress;anisotropic;triangular;tri mesh";

            // ── Input ports ───────────────────────────────────────────────────────────────

            // Port index 0: the solved FEA model. Must contain only triangular shell elements.
            InputParameterManager.AddParameter<IModel>("Input model",
                "Solved FEA model with triangular shell elements whose stress results will guide the mesh deformation.",
                ParameterAccess.Item);

            // Port index 1: number of relaxation passes.
            // Triangular meshes typically need more iterations than quad meshes
            // because each step moves nodes a smaller fraction of the way (due to Lambda damping).
            InputParameterManager.AddParameter<SyneraInt>("Iterations",
                "Number of relaxation iterations. Higher values produce smoother, more stress-aligned meshes. " +
                "Triangular meshes typically need more iterations than quad meshes to reach a stable result.",
                ParameterAccess.Item);

            // Port index 2: how strongly von Mises stress pulls nodes toward stressed regions.
            // Alpha = 0 → uniform smoothing. Higher values → denser mesh near stress peaks.
            InputParameterManager.AddParameter<SyneraInt>("Alpha",
                "Controls how strongly von Mises stress influences node movement relative to the geometric center. " +
                "Alpha = 0 gives uniform Laplacian smoothing; higher values pull nodes toward high-stress regions.",
                ParameterAccess.Item);

            // Port index 3: how strongly principal stress directions steer node movement.
            // Beta = 0 → isotropic smoothing (same as FieldCreation).
            // Higher Beta → elements elongate along principal stress trajectories.
            // Step size is always preserved — only direction changes — so inversion is not possible.
            InputParameterManager.AddParameter<SyneraDouble>("Beta",
                "Controls how strongly principal stress directions influence node movement. " +
                "Beta = 0 disables directional anisotropy; higher values rotate the displacement " +
                "progressively toward the principal stress axes without changing the step size, " +
                "so mesh inversion is not possible.",
                ParameterAccess.Item);

            // Port index 4: the damping factor — how far each node moves per iteration.
            // Lambda = 1.0 → full step (node moves all the way to the target position).
            // Lambda = 0.5 → half step (node moves halfway, then recalculates next iteration).
            // Lower Lambda → slower but smoother convergence, safer for triangular meshes
            //   that are prone to becoming slivers when pushed too hard in one direction.
            // Default is 1.0 (no damping). Reduce to 0.5 if the result looks distorted.
            InputParameterManager.AddParameter<SyneraDouble>("Lambda",
                "Damping factor controlling what fraction of the computed displacement is applied each iteration. " +
                "Lambda = 1.0 applies the full step (no damping); Lambda = 0.5 applies half the step. " +
                "Lower values produce smoother convergence and prevent element distortion in triangular meshes. " +
                "Valid range: (0, 1]. Defaults to 1.0.",
                ParameterAccess.Item,
                new SyneraDouble(1.0)); // default value — port can be left unconnected

            // ── Output ports ──────────────────────────────────────────────────────────────

            // Port index 0: major principal stress direction vector per element.
            OutputParameterManager.AddParameter<Vector3D>("Max Principle Stress",
                "Maximum principal stress vectors per element.",
                ParameterAccess.List);

            // Port index 1: minor principal stress direction vector per element.
            OutputParameterManager.AddParameter<Vector3D>("Min Principle Stress",
                "Minimum principal stress vectors per element.",
                ParameterAccess.List);

            // Port index 2: von Mises stress averaged per node (from surrounding elements).
            OutputParameterManager.AddParameter<SyneraDouble>("Von-Mises Stress",
                "Von Mises stress value averaged per node.",
                ParameterAccess.List);

            // Port index 3: original node positions before any deformation.
            OutputParameterManager.AddParameter<Point3D>("Vertices",
                "Original node locations before deformation.",
                ParameterAccess.List);

            // Port index 4: node positions after the stress-guided relaxation.
            OutputParameterManager.AddParameter<Point3D>("Deformed Vertices",
                "Node locations after stress-guided anisotropic deformation.",
                ParameterAccess.List);

            // Port index 5: the rebuilt mesh using the deformed vertex positions.
            OutputParameterManager.AddParameter<IMesh>("Mesh",
                "Resulting deformed mesh aligned with the principal stress field.",
                ParameterAccess.Item);
        }

        /// <summary>
        /// Converts a 6-component stress tensor into a single von Mises stress scalar.
        ///
        /// Von Mises stress combines normal and shear stress components into one value
        /// that represents the overall stress intensity at a point. It is commonly used
        /// to predict where a structure is likely to yield or fail.
        ///
        /// Input order: [Sxx, Syy, Szz, Sxy, Syz, Sxz]
        /// </summary>
        private static double CalculateVonMisesStress(double[] stressTensor)
        {
            // Extract the six independent components from the symmetric stress tensor.
            double sxx = stressTensor[0], syy = stressTensor[1], szz = stressTensor[2];
            double sxy = stressTensor[3], syz = stressTensor[4], sxz = stressTensor[5];

            // Normal stress contribution: measures how unevenly the material is being stretched.
            double normalPart = (sxx - syy) * (sxx - syy)
                              + (syy - szz) * (syy - szz)
                              + (szz - sxx) * (szz - sxx);

            // Shear stress contribution: measures how much the material is being twisted.
            double shearPart = 6.0 * (sxy * sxy + syz * syz + sxz * sxz);

            return Math.Sqrt(0.5 * (normalPart + shearPart));
        }

        /// <summary>
        /// Moves interior triangular mesh nodes iteratively to align with the stress field.
        ///
        /// WHY JACOBI FOR TRIANGULAR MESHES:
        ///   A standard Gauss-Seidel update applies each node's new position immediately,
        ///   so later nodes in the loop see a mix of old and new neighbor positions.
        ///   For quad meshes (~4 neighbors at 90°), this is mild and stable.
        ///   For triangular meshes (~6 neighbors at 60°), the denser connectivity amplifies
        ///   this asymmetry, causing chaotic-looking deformation after several iterations.
        ///
        ///   The Jacobi approach fixes this by:
        ///     1. Computing ALL new positions using only the positions from the START of the iteration.
        ///     2. Applying ALL new positions at the END of the iteration simultaneously.
        ///   This means every node sees the same consistent snapshot of its neighbors.
        ///
        /// WHY DAMPING (LAMBDA):
        ///   Triangles resist elongation geometrically — pushing them too far in one direction
        ///   per iteration creates sliver elements. Lambda limits each step to a fraction of the
        ///   computed displacement, allowing gradual convergence without element degradation.
        ///
        /// ALGORITHM PER NODE PER ITERATION:
        ///   Stage 1 — Von Mises weighted centroid (read from frozen snapshot):
        ///     Compute weighted average of neighbor positions. Neighbors with higher stress
        ///     get more weight (controlled by Alpha).
        ///   Stage 2 — Decompose displacement into stress-direction components:
        ///     Split the node→centroid vector into major stress, minor stress, and residual parts.
        ///   Stage 3 — Amplify stress-aligned components, preserve step length:
        ///     Scale up the stress-direction parts by Beta, then rescale the whole vector
        ///     back to the original length. Direction changes; magnitude does not.
        ///   Stage 4 — Apply damping:
        ///     Move only Lambda fraction of the computed displacement.
        /// </summary>
        private static void DeformNodes(
            List<INode> allNodes,           // All nodes in the mesh
            int iterations,                 // Number of relaxation passes
            double alpha,                   // Von Mises stress influence strength
            double beta,                    // Principal stress direction influence strength
            double lambda,                  // Damping factor — fraction of displacement applied per iteration
            Dictionary<int, double> nodeVonMisesStress,             // Von Mises stress value per node index
            Dictionary<int, List<int>> nodeToNeighborsMap,          // Neighbor node indices per node
            HashSet<int> boundaryNodeIndices,                        // Nodes locked in place
            Dictionary<int, (Vector3D majorDirection, Vector3D minorDirection)> nodePrincipalDirections) // Stress directions per node
        {
            // Normalize von Mises stress to [0, 1] so Alpha is scale-independent.
            double minStress = nodeVonMisesStress.Values.Min();
            double maxStress = nodeVonMisesStress.Values.Max();
            double stressRange = maxStress - minStress + 1e-6; // epsilon prevents division by zero

            // The Jacobi read buffer: stores node positions at the START of each iteration.
            // All centroid computations within one iteration read from here, never from
            // partially-updated positions. This is the key difference from Gauss-Seidel.
            Dictionary<int, Point3D> currentPositions = allNodes.ToDictionary(n => n.Index, n => n.Location);

            for (int iteration = 0; iteration < iterations; iteration++)
            {
                // The Jacobi write buffer: all new positions are stored here during the pass
                // and applied to the actual nodes only after the full pass is complete.
                Dictionary<int, Point3D> newPositions = new Dictionary<int, Point3D>(allNodes.Count);

                foreach (INode node in allNodes)
                {
                    // Boundary nodes stay fixed — copy their current position unchanged.
                    if (boundaryNodeIndices.Contains(node.Index))
                    {
                        newPositions[node.Index] = currentPositions[node.Index];
                        continue;
                    }

                    List<int> neighbors = nodeToNeighborsMap[node.Index];
                    if (neighbors.Count == 0)
                    {
                        newPositions[node.Index] = currentPositions[node.Index];
                        continue;
                    }

                    // Retrieve the averaged principal stress directions for this node.
                    // These are computed by averaging over all elements connected to the node.
                    nodePrincipalDirections.TryGetValue(node.Index, out var principalDirs);
                    Vector3D majorDir = principalDirs.majorDirection;
                    Vector3D minorDir = principalDirs.minorDirection;

                    double majorMagnitude = majorDir.Length;
                    double minorMagnitude = minorDir.Length;

                    // Per-node local normalization: the dominant stress direction always scores 1.0.
                    // Using a global maximum instead would cause stress peaks to flatten
                    // the contribution of all other nodes to near zero.
                    double localMax = Math.Max(majorMagnitude, minorMagnitude) + 1e-6;
                    double normalizedMajorMagnitude = majorMagnitude / localMax;
                    double normalizedMinorMagnitude = minorMagnitude / localMax;

                    // Unit vectors pointing along each principal stress direction.
                    // Fall back to world axes for nodes with negligible stress (e.g. free corners).
                    Vector3D majorUnit = majorMagnitude > 1e-10 ? majorDir.Normalized() : new Vector3D(1, 0, 0);
                    Vector3D minorUnit = minorMagnitude > 1e-10 ? minorDir.Normalized() : new Vector3D(0, 1, 0);

                    // ── Stage 1: Von Mises weighted centroid (Jacobi read) ────────────────────
                    // Read neighbor positions from currentPositions — the frozen snapshot.
                    // This ensures all nodes in this iteration use positions from the same moment.
                    Point3D weightedPositionSum = new Point3D(0, 0, 0);
                    double totalWeight = 0.0;

                    foreach (int neighborIdx in neighbors)
                    {
                        // The weight formula gives higher influence to high-stress neighbors.
                        // At Alpha = 0, all neighbors have equal weight = 1.0 (pure Laplacian).
                        // At Alpha = 5, a max-stress neighbor has weight = 6.0, min-stress = 1.0.
                        double normalizedStress = (nodeVonMisesStress[neighborIdx] - minStress) / stressRange;
                        double vonMisesWeight = 1.0 + alpha * normalizedStress;

                        weightedPositionSum += vonMisesWeight * currentPositions[neighborIdx];
                        totalWeight += vonMisesWeight;
                    }

                    if (totalWeight <= 1e-10)
                    {
                        newPositions[node.Index] = currentPositions[node.Index];
                        continue;
                    }

                    // This is where the node would end up under pure von Mises weighted smoothing.
                    Point3D centroid = weightedPositionSum / totalWeight;

                    // ── Stage 2: Decompose displacement into principal stress components ─────────
                    // The displacement vector points from the node's current position to the centroid.
                    Point3D oldPosition = currentPositions[node.Index];
                    Vector3D delta = Point3D.Subtract(centroid, oldPosition);

                    // Project delta onto each principal stress direction.
                    // The scalar dot product gives the signed length of the projection.
                    double dMajor = delta * majorUnit;
                    double dMinor = delta * minorUnit;

                    // The residual is the part of delta that doesn't align with either stress direction.
                    // For planar shell meshes this is mostly the through-thickness component.
                    Vector3D deltaRest = delta - dMajor * majorUnit - dMinor * minorUnit;

                    // ── Stage 3: Steer direction, preserve step length ─────────────────────────
                    // Amplify the stress-aligned components by Beta.
                    // At Beta = 0, anisotropicDirection = delta (no steering).
                    // At Beta > 0, the major/minor components grow relative to deltaRest,
                    //   so the resulting direction tilts toward the principal stress axes.
                    Vector3D anisotropicDirection = deltaRest
                        + (1.0 + beta * normalizedMajorMagnitude) * dMajor * majorUnit
                        + (1.0 + beta * normalizedMinorMagnitude) * dMinor * minorUnit;

                    // Rescale back to the original delta length so Beta only changes direction.
                    // Without this, high Beta would move nodes further than their centroid,
                    // causing overshoot and possible mesh inversion.
                    double deltaLength = delta.Length;
                    double anisotropicLength = anisotropicDirection.Length;

                    Vector3D anisotropicDelta;
                    if (anisotropicLength > 1e-10 && deltaLength > 1e-10)
                        anisotropicDelta = (deltaLength / anisotropicLength) * anisotropicDirection;
                    else
                        anisotropicDelta = delta; // fallback: use the unsteered displacement

                    // ── Stage 4: Apply with damping ──────────────────────────────────────────
                    // Lambda < 1.0 means the node moves only part of the way to its target.
                    // This prevents triangular elements from becoming thin slivers by limiting
                    // how much the geometry changes in a single iteration.
                    newPositions[node.Index] = oldPosition + lambda * anisotropicDelta;
                }

                // ── Jacobi write: apply all new positions at once ─────────────────────────────
                // Only after ALL nodes have been computed do we update their actual positions.
                // This guarantees that no node in this iteration influenced another node's result.
                foreach (INode node in allNodes)
                {
                    if (newPositions.TryGetValue(node.Index, out Point3D newPos))
                    {
                        node.Location = newPos;
                        // Also update the read buffer so the next iteration starts from
                        // the positions we just applied.
                        currentPositions[node.Index] = newPos;
                    }
                }
            }
        }

        /// <summary>
        /// Called by Synera whenever the node needs to recompute.
        /// Reads and validates inputs, runs the relaxation, and writes outputs.
        /// </summary>
        protected override void SolveInstance(IDataAccess dataAccess)
        {
            // ── Read inputs from the connected ports ──────────────────────────────────────
            bool isDataSuccess = dataAccess.GetData(0, out IModel model);
            isDataSuccess &= dataAccess.GetData(1, out SyneraInt iterations);
            isDataSuccess &= dataAccess.GetData(2, out SyneraInt alpha);
            isDataSuccess &= dataAccess.GetData(3, out double beta);
            isDataSuccess &= dataAccess.GetData(4, out double lambda);

            // If any required port is disconnected or has no data, abort.
            if (!isDataSuccess)
                return;

            // ── Validate inputs ───────────────────────────────────────────────────────────
            // Collect all validation errors before returning so the user sees all issues at once.
            bool isValid = true;

            if (!model.IsValid)
            {
                AddError(0, "Invalid model.");
                isValid = false;
            }
            else if (!model.HasShellElementsOnly())
            {
                AddError(0, "A model containing shell elements only is expected.");
                isValid = false;
            }
            else
            {
                // Check that every shell element is triangular (3 corner nodes).
                // This node is specifically designed for triangular topology.
                // Quad meshes should use QuadMeshFD instead.
                IEnumerable<IShellElement> shellElements = model.Elements.OfType<IShellElement>();
                bool hasTrianglesOnly = shellElements.All(elem =>
                    elem.NodeIndices.Except(elem.GetMidNodeIndices()).Count() == 3);

                if (!hasTrianglesOnly)
                {
                    AddError(0, "This node expects a triangular shell mesh. Use Quad Mesh FD for quad elements.");
                    isValid = false;
                }
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

            // Lambda must be strictly positive (0 means no movement — useless)
            // and at most 1.0 (full step — values above 1 would overshoot the target).
            if (lambda <= 0 || lambda > 1.0)
            {
                AddError(4, "Lambda must be in the range (0, 1]. Use 0.5 as a safe starting value.");
                isValid = false;
            }

            if (!isValid)
                return;

            // ── Step 1: Compute principal stress vectors per element ──────────────────────
            // StressStreamlineBuilder performs eigenvalue decomposition on the stress tensor
            // of each shell element to find its two in-plane principal stress directions.
            Progress progress = this.CreateProgress();
            StressStreamlineBuilder streamlineBuilder = new StressStreamlineBuilder(model);

            Progress stressProgress = progress.CreateSubtask(0.5);
            Dictionary<int, (Vector3D majorVector, Vector3D minorVector)> elementPrincipalStress =
                streamlineBuilder.GetPrincipalStressValues(stressProgress);

            // Split into separate lists for the output ports.
            List<Vector3D> majorStressVectors = elementPrincipalStress.Values.Select(v => v.majorVector).ToList();
            List<Vector3D> minorStressVectors  = elementPrincipalStress.Values.Select(v => v.minorVector).ToList();

            // ── Step 2: Build mesh topology lookup maps ───────────────────────────────────
            IEnumerable<IShellElement> elements = model.Elements.OfType<IShellElement>();
            List<INode> nodes = model.Nodes.ToList();

            // Map: node index → indices of elements that share this node.
            // Only corner nodes are included — mid-side nodes of higher-order elements are excluded.
            Dictionary<int, List<int>> nodeToElementsMap = new Dictionary<int, List<int>>();
            foreach (IShellElement element in elements)
            {
                List<int> cornerNodeIndices = element.NodeIndices.Except(element.GetMidNodeIndices()).ToList();
                foreach (int nodeIndex in cornerNodeIndices)
                {
                    if (!nodeToElementsMap.ContainsKey(nodeIndex))
                        nodeToElementsMap[nodeIndex] = new List<int>();
                    nodeToElementsMap[nodeIndex].Add(element.Index);
                }
            }

            // Map: node index → indices of directly connected neighbor nodes (edge-adjacent).
            // Two nodes are neighbors if they share a mesh edge.
            Dictionary<int, List<int>> nodeToNeighborsMap = new Dictionary<int, List<int>>();
            foreach (IndexPair edge in model.Edges)
            {
                int a = edge.A;
                int b = edge.B;

                // Skip degenerate zero-length edges.
                if (a == b)
                    continue;

                if (!nodeToNeighborsMap.ContainsKey(a)) nodeToNeighborsMap[a] = new List<int>();
                if (!nodeToNeighborsMap.ContainsKey(b)) nodeToNeighborsMap[b] = new List<int>();

                // Each edge connects two nodes — add each as the other's neighbor (bidirectional).
                if (!nodeToNeighborsMap[a].Contains(b)) nodeToNeighborsMap[a].Add(b);
                if (!nodeToNeighborsMap[b].Contains(a)) nodeToNeighborsMap[b].Add(a);
            }

            // ── Step 3: Compute per-node von Mises stress ─────────────────────────────────
            // FEA provides stress per element. We average it onto each node
            // by combining the stress values from all elements connected to that node.
            IResultsAtTimeStep elementStressResults = model.Results
                .First(r => r.Name.Equals(ResultType.ElementStressTensor.Name))
                .Values.First()  // first load case
                .Values.First(); // first time step

            Dictionary<int, double> nodeVonMisesStress = new Dictionary<int, double>(nodes.Count);
            foreach (INode node in nodes)
            {
                if (nodeToElementsMap.TryGetValue(node.Index, out List<int> connectedElements))
                {
                    double stressSum = connectedElements
                        .Sum(elemIdx => CalculateVonMisesStress(elementStressResults[elemIdx]));
                    nodeVonMisesStress[node.Index] = stressSum / connectedElements.Count;
                }
                else
                {
                    // Nodes not connected to any element carry zero stress.
                    nodeVonMisesStress[node.Index] = 0.0;
                }
            }

            // ── Step 4: Average principal stress directions per node ───────────────────────
            // Element-level principal stress vectors are averaged over all elements
            // connected to each node, giving a smooth directional field at node locations.
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
                    double inv = 1.0 / contributingElements;
                    nodePrincipalDirections[node.Index] = (inv * majorSum, inv * minorSum);
                }
            }

            // ── Step 5: Identify boundary nodes that must stay fixed ──────────────────────
            // Nodes on naked (open boundary) edges define the mesh outline.
            // Keeping them fixed ensures the deformation does not shrink or distort the footprint.
            IEnumerable<int> nakedEdgeIndices = streamlineBuilder.Mesh.GetNakedEdgeIndices();
            HashSet<int> boundaryNodeIndices = new HashSet<int>();
            foreach (int edgeIdx in nakedEdgeIndices)
            {
                IndexPair edge = model.GetNodeIndicesForEdge(edgeIdx);
                boundaryNodeIndices.AddRange(edge.ToList());
            }

            // Record original positions before any movement for the output port.
            List<Point3D> originalVertices = nodes.Select(n => n.Location).ToList();

            // ── Step 6: Run the Jacobi-based stress-guided relaxation ─────────────────────
            // This modifies node positions in-place using the Jacobi update scheme.
            DeformNodes(
                nodes,
                iterations,
                alpha,
                beta,
                lambda,
                nodeVonMisesStress,
                nodeToNeighborsMap,
                boundaryNodeIndices,
                nodePrincipalDirections);

            // Collect updated node positions after all iterations are complete.
            List<Point3D> deformedVertices = nodes.Select(n => n.Location).ToList();

            // ── Step 7: Rebuild the mesh using the deformed vertex positions ───────────────
            // Mesh connectivity (which nodes form which faces) is unchanged.
            // Only the spatial positions of the nodes are different.
            List<MeshFace> meshFaces = elements
                .Select(elem => new MeshFace(elem.NodeIndices.Select(i => i - 1).ToArray()))
                .ToList();

            IMeshKernel meshKernel = Application.Current.KernelManager.Get<IMeshKernel>();
            IMesh deformedMesh = meshKernel.CreateFromVerticesAndFaces(deformedVertices, meshFaces);

            // ── Step 8: Write all outputs to the connected ports ──────────────────────────
            dataAccess.SetListData(0, majorStressVectors);
            dataAccess.SetListData(1, minorStressVectors);
            dataAccess.SetListData(2, nodeVonMisesStress.Values.ToList());
            dataAccess.SetListData(3, originalVertices);
            dataAccess.SetListData(4, deformedVertices);
            dataAccess.SetData(5, deformedMesh);
        }
    }
}
