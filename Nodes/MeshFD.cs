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
    [Guid("7B2C4D5E-8F9A-6B7C-8D9E-0F1A2B3C4D5E")]
    public sealed class MeshFD : Node
    {
        public MeshFD() : base("Mesh FD")
        {
            Category = WobisCategories.Aida;
            Subcategory = WobisSubcategories.Field;
            Description = "Deforms a shell mesh (quad or triangular) guided by both von Mises stress " +
                          "magnitude and principal stress directions. Mesh type is detected automatically: " +
                          "quad meshes use a fast Gauss-Seidel update; triangular meshes use a synchronous " +
                          "Jacobi update with per-iteration damping (Lambda) to prevent element distortion.";
            GuiPriority = 2;
            Keywords = "deform mesh;stress-guided;principal stress;anisotropic;quad;triangular";

            InputParameterManager.AddParameter<IModel>("Input model",
                "Solved FEA model with shell elements (quad or triangular) whose stress results will guide the deformation.",
                ParameterAccess.Item);

            InputParameterManager.AddParameter<SyneraInt>("Iterations",
                "Number of relaxation iterations. Higher values produce smoother, more stress-aligned meshes.",
                ParameterAccess.Item);

            InputParameterManager.AddParameter<SyneraInt>("Alpha",
                "Controls how strongly von Mises stress influences node movement. " +
                "Alpha = 0 gives uniform Laplacian smoothing; higher values pull nodes toward high-stress regions.",
                ParameterAccess.Item);

            InputParameterManager.AddParameter<SyneraDouble>("Beta",
                "Controls how strongly principal stress directions influence node movement. " +
                "Beta = 0 disables directional anisotropy; higher values rotate the displacement " +
                "toward the principal stress axes without changing the step size.",
                ParameterAccess.Item);

            InputParameterManager.AddParameter<SyneraDouble>("Lambda",
                "Damping factor for triangular meshes — controls what fraction of the computed " +
                "displacement is applied each iteration. Valid range: (0, 1]. " +
                "Lambda = 1.0 applies the full step; Lambda = 0.5 applies half the step. " +
                "Ignored for quad meshes (always treated as 1.0). Defaults to 1.0.",
                ParameterAccess.Item,
                new SyneraDouble(1.0));

            OutputParameterManager.AddParameter<Vector3D>("Max Principle Stress",
                "Maximum principal stress vectors per element.",
                ParameterAccess.List);

            OutputParameterManager.AddParameter<Vector3D>("Min Principle Stress",
                "Minimum principal stress vectors per element.",
                ParameterAccess.List);

            OutputParameterManager.AddParameter<SyneraDouble>("Von-Mises Stress",
                "Von Mises stress value averaged per node.",
                ParameterAccess.List);

            OutputParameterManager.AddParameter<Point3D>("Vertices",
                "Original node locations before deformation.",
                ParameterAccess.List);

            OutputParameterManager.AddParameter<Point3D>("Deformed Vertices",
                "Node locations after stress-guided anisotropic deformation.",
                ParameterAccess.List);

            OutputParameterManager.AddParameter<IMesh>("Mesh",
                "Resulting deformed mesh aligned with the principal stress field.",
                ParameterAccess.Item);
        }

        /// <summary>
        /// Computes the von Mises stress scalar from a 6-component stress tensor.
        /// Components are expected in the order [Sxx, Syy, Szz, Sxy, Syz, Sxz].
        /// </summary>
        private static double CalculateVonMisesStress(double[] stressTensor)
        {
            double sxx = stressTensor[0], syy = stressTensor[1], szz = stressTensor[2];
            double sxy = stressTensor[3], syz = stressTensor[4], sxz = stressTensor[5];

            double normalPart = (sxx - syy) * (sxx - syy)
                              + (syy - szz) * (syy - szz)
                              + (szz - sxx) * (szz - sxx);

            double shearPart = 6.0 * (sxy * sxy + syz * syz + sxz * sxz);

            return Math.Sqrt(0.5 * (normalPart + shearPart));
        }

        /// <summary>
        /// Unified anisotropic stress-guided relaxation for both quad and triangular meshes.
        ///
        /// UPDATE SCHEME — controlled by useJacobi:
        ///
        ///   Gauss-Seidel (useJacobi = false) — quad meshes:
        ///     Node positions are updated in-place immediately within the iteration loop.
        ///     Subsequent nodes read the already-updated positions of earlier neighbors.
        ///     Converges faster per iteration. Stable for quads (~4 neighbors at 90°).
        ///
        ///   Jacobi (useJacobi = true) — triangular meshes:
        ///     All new positions are computed from a frozen snapshot taken at the start of
        ///     each iteration and applied simultaneously at the end. No node sees a
        ///     partially-updated neighborhood mid-iteration. Prevents the chaotic propagation
        ///     that occurs in dense triangular connectivity (~6 neighbors at 60°).
        ///
        /// DAMPING — controlled by lambda:
        ///   Each node moves only a fraction (lambda) of the computed anisotropic displacement.
        ///   lambda = 1.0 → full step (quads). lambda < 1.0 → partial step (triangles),
        ///   preventing sliver formation by limiting per-iteration geometric change.
        ///
        /// ALGORITHM PER NODE PER ITERATION:
        ///   1. Compute a von Mises weighted Laplacian centroid from neighbor positions
        ///      (Gauss-Seidel: live positions; Jacobi: frozen snapshot positions).
        ///   2. Decompose the displacement vector (node → centroid) into components along
        ///      the major and minor principal stress directions and a residual.
        ///   3. Amplify the stress-aligned components by Beta, then rescale the entire
        ///      direction vector back to the original displacement length — only direction
        ///      changes, step size is never increased, preventing overshoot.
        ///   4. Apply Lambda fraction of the resulting displacement to the node position.
        ///   5. (Jacobi only) Write all new positions simultaneously after the full pass.
        /// </summary>
        private static void DeformNodes(
            List<INode> allNodes,
            int iterations,
            double alpha,
            double beta,
            double lambda,
            bool useJacobi,
            Dictionary<int, double> nodeVonMisesStress,
            Dictionary<int, List<int>> nodeToNeighborsMap,
            HashSet<int> boundaryNodeIndices,
            Dictionary<int, (Vector3D majorDirection, Vector3D minorDirection)> nodePrincipalDirections)
        {
            double minStress = nodeVonMisesStress.Values.Min();
            double maxStress = nodeVonMisesStress.Values.Max();
            double stressRange = maxStress - minStress + 1e-6;

            // nodeDict: live lookup for Gauss-Seidel (reflects in-place updates immediately).
            Dictionary<int, INode> nodeDict = allNodes.ToDictionary(n => n.Index);

            // currentPositions: frozen snapshot for Jacobi (never updated mid-iteration).
            Dictionary<int, Point3D> currentPositions = allNodes.ToDictionary(n => n.Index, n => n.Location);

            for (int iteration = 0; iteration < iterations; iteration++)
            {
                // Jacobi write buffer — allocated only when needed.
                Dictionary<int, Point3D> newPositions = useJacobi
                    ? new Dictionary<int, Point3D>(allNodes.Count)
                    : null;

                foreach (INode node in allNodes)
                {
                    // Boundary nodes remain fixed — copy position unchanged for Jacobi.
                    if (boundaryNodeIndices.Contains(node.Index))
                    {
                        if (useJacobi) newPositions[node.Index] = currentPositions[node.Index];
                        continue;
                    }

                    List<int> neighbors = nodeToNeighborsMap[node.Index];
                    if (neighbors.Count == 0)
                    {
                        if (useJacobi) newPositions[node.Index] = currentPositions[node.Index];
                        continue;
                    }

                    // Retrieve per-node averaged principal stress directions.
                    nodePrincipalDirections.TryGetValue(node.Index, out var principalDirs);
                    Vector3D majorDir = principalDirs.majorDirection;
                    Vector3D minorDir = principalDirs.minorDirection;

                    double majorMagnitude = majorDir.Length;
                    double minorMagnitude = minorDir.Length;

                    // Per-node local normalization: dominant direction always scores 1.0,
                    // subordinate direction scores proportionally. Prevents global stress
                    // peaks from collapsing the directional contribution of all other nodes.
                    double localMax = Math.Max(majorMagnitude, minorMagnitude) + 1e-6;
                    double normalizedMajorMagnitude = majorMagnitude / localMax;
                    double normalizedMinorMagnitude = minorMagnitude / localMax;

                    // Unit vectors along principal stress directions.
                    // Fall back to world axes if magnitudes are negligible (stress-free node).
                    Vector3D majorUnit = majorMagnitude > 1e-10 ? majorDir.Normalized() : new Vector3D(1, 0, 0);
                    Vector3D minorUnit = minorMagnitude > 1e-10 ? minorDir.Normalized() : new Vector3D(0, 1, 0);

                    // ── Step 1: Von Mises weighted centroid ───────────────────────────────────
                    // Gauss-Seidel reads live positions from nodeDict (some already updated).
                    // Jacobi reads from the frozen snapshot in currentPositions (all from start of iteration).
                    Point3D weightedPositionSum = new Point3D(0, 0, 0);
                    double totalWeight = 0.0;

                    foreach (int neighborIdx in neighbors)
                    {
                        Point3D neighborPos = useJacobi
                            ? currentPositions[neighborIdx]
                            : nodeDict[neighborIdx].Location;

                        double normalizedStress = (nodeVonMisesStress[neighborIdx] - minStress) / stressRange;
                        double vonMisesWeight = 1.0 + alpha * normalizedStress;

                        weightedPositionSum += vonMisesWeight * neighborPos;
                        totalWeight += vonMisesWeight;
                    }

                    if (totalWeight <= 1e-10)
                    {
                        if (useJacobi) newPositions[node.Index] = currentPositions[node.Index];
                        continue;
                    }

                    Point3D centroid = weightedPositionSum / totalWeight;

                    // ── Step 2: Decompose displacement into principal stress components ─────────
                    Point3D oldPosition = useJacobi ? currentPositions[node.Index] : node.Location;
                    Vector3D delta = Point3D.Subtract(centroid, oldPosition);

                    double dMajor = delta * majorUnit;
                    double dMinor = delta * minorUnit;
                    Vector3D deltaRest = delta - dMajor * majorUnit - dMinor * minorUnit;

                    // ── Step 3: Steer direction toward principal stress axes ───────────────────
                    // Beta amplifies stress-aligned components. The result is then rescaled to
                    // the original delta length so the step size is never increased — Beta
                    // rotates the direction only, overshoot and inversion are not possible.
                    Vector3D anisotropicDirection = deltaRest
                        + (1.0 + beta * normalizedMajorMagnitude) * dMajor * majorUnit
                        + (1.0 + beta * normalizedMinorMagnitude) * dMinor * minorUnit;

                    double deltaLength = delta.Length;
                    double anisotropicLength = anisotropicDirection.Length;

                    Vector3D anisotropicDelta;
                    if (anisotropicLength > 1e-10 && deltaLength > 1e-10)
                        anisotropicDelta = (deltaLength / anisotropicLength) * anisotropicDirection;
                    else
                        anisotropicDelta = delta;

                    // ── Step 4: Apply with damping ────────────────────────────────────────────
                    // Lambda blends between no movement (0) and the full anisotropic step (1).
                    Point3D newPosition = oldPosition + lambda * anisotropicDelta;

                    if (useJacobi)
                        newPositions[node.Index] = newPosition;
                    else
                        node.Location = newPosition;
                }

                // ── Jacobi write: apply all new positions simultaneously ──────────────────────
                // This is the key difference from Gauss-Seidel — no node sees another node's
                // updated position until the next iteration begins.
                if (useJacobi)
                {
                    foreach (INode node in allNodes)
                    {
                        if (newPositions.TryGetValue(node.Index, out Point3D newPos))
                        {
                            node.Location = newPos;
                            currentPositions[node.Index] = newPos;
                        }
                    }
                }
            }
        }

        protected override void SolveInstance(IDataAccess dataAccess)
        {
            // ── Read inputs ───────────────────────────────────────────────────────────────
            bool isDataSuccess = dataAccess.GetData(0, out IModel model);
            isDataSuccess &= dataAccess.GetData(1, out SyneraInt iterations);
            isDataSuccess &= dataAccess.GetData(2, out SyneraInt alpha);
            isDataSuccess &= dataAccess.GetData(3, out double beta);
            isDataSuccess &= dataAccess.GetData(4, out double lambda);

            if (!isDataSuccess)
                return;

            // ── Validate inputs ───────────────────────────────────────────────────────────
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

            if (lambda <= 0 || lambda > 1.0)
            {
                AddError(4, "Lambda must be in the range (0, 1]. Use 0.5 as a safe starting value for triangular meshes.");
                isValid = false;
            }

            if (!isValid)
                return;

            // ── Detect mesh type and configure update strategy ────────────────────────────
            IEnumerable<IShellElement> shellElements = model.Elements.OfType<IShellElement>();

            // A mesh is purely triangular only if every element has exactly 3 corner nodes.
            // Any other case (all quads, or a mix of quads and triangles) is treated as
            // quad-dominant and processed with the faster Gauss-Seidel update.
            bool isTriangular = shellElements.All(elem =>
                elem.NodeIndices.Except(elem.GetMidNodeIndices()).Count() == 3);

            bool isMixed = !isTriangular && shellElements.Any(elem =>
                elem.NodeIndices.Except(elem.GetMidNodeIndices()).Count() == 3);

            // Pure tri  → Jacobi update + Lambda damping (stable for dense tri connectivity).
            // Quad / Mixed → Gauss-Seidel update, Lambda forced to 1.0 (fast, stable for quads).
            bool useJacobi = isTriangular;
            double effectiveLambda = isTriangular ? lambda : 1.0;

            if (isTriangular)
                AddNotification("Triangular mesh detected — Jacobi update with Lambda damping applied.");
            else if (isMixed)
                AddNotification("Mixed quad/triangular mesh detected — treated as quad-dominant. Gauss-Seidel update applied. Lambda ignored.");
            else
                AddNotification("Quad mesh detected — Gauss-Seidel update applied. Lambda ignored.");

            // ── Compute principal stress vectors ──────────────────────────────────────────
            Progress progress = this.CreateProgress();
            StressStreamlineBuilder streamlineBuilder = new StressStreamlineBuilder(model);

            Progress stressProgress = progress.CreateSubtask(0.5);
            Dictionary<int, (Vector3D majorVector, Vector3D minorVector)> elementPrincipalStress =
                streamlineBuilder.GetPrincipalStressValues(stressProgress);

            List<Vector3D> majorStressVectors = elementPrincipalStress.Values.Select(v => v.majorVector).ToList();
            List<Vector3D> minorStressVectors  = elementPrincipalStress.Values.Select(v => v.minorVector).ToList();

            // ── Build mesh topology maps ──────────────────────────────────────────────────
            List<INode> nodes = model.Nodes.ToList();

            // Map: node index → indices of all elements that share this node.
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

            // Map: node index → indices of all directly connected (edge-adjacent) nodes.
            Dictionary<int, List<int>> nodeToNeighborsMap = new Dictionary<int, List<int>>();
            foreach (IndexPair edge in model.Edges)
            {
                int a = edge.A;
                int b = edge.B;

                if (a == b)
                    continue;

                if (!nodeToNeighborsMap.ContainsKey(a)) nodeToNeighborsMap[a] = new List<int>();
                if (!nodeToNeighborsMap.ContainsKey(b)) nodeToNeighborsMap[b] = new List<int>();

                if (!nodeToNeighborsMap[a].Contains(b)) nodeToNeighborsMap[a].Add(b);
                if (!nodeToNeighborsMap[b].Contains(a)) nodeToNeighborsMap[b].Add(a);
            }

            // ── Compute per-node von Mises stress ─────────────────────────────────────────
            IResultsAtTimeStep elementStressResults = model.Results
                .First(r => r.Name.Equals(ResultType.ElementStressTensor.Name))
                .Values.First()
                .Values.First();

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
                    nodeVonMisesStress[node.Index] = 0.0;
                }
            }

            // ── Average principal stress directions per node ───────────────────────────────
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

            // ── Identify fixed boundary nodes ─────────────────────────────────────────────
            IEnumerable<int> nakedEdgeIndices = streamlineBuilder.Mesh.GetNakedEdgeIndices();
            HashSet<int> boundaryNodeIndices = new HashSet<int>();
            foreach (int edgeIdx in nakedEdgeIndices)
            {
                IndexPair edge = model.GetNodeIndicesForEdge(edgeIdx);
                boundaryNodeIndices.AddRange(edge.ToList());
            }

            List<Point3D> originalVertices = nodes.Select(n => n.Location).ToList();

            // ── Run unified stress-guided relaxation ──────────────────────────────────────
            DeformNodes(
                nodes,
                iterations,
                alpha,
                beta,
                effectiveLambda,
                useJacobi,
                nodeVonMisesStress,
                nodeToNeighborsMap,
                boundaryNodeIndices,
                nodePrincipalDirections);

            List<Point3D> deformedVertices = nodes.Select(n => n.Location).ToList();

            // ── Build the output mesh ─────────────────────────────────────────────────────
            List<MeshFace> meshFaces = shellElements
                .Select(elem => new MeshFace(elem.NodeIndices.Select(i => i - 1).ToArray()))
                .ToList();

            IMeshKernel meshKernel = Application.Current.KernelManager.Get<IMeshKernel>();
            IMesh deformedMesh = meshKernel.CreateFromVerticesAndFaces(deformedVertices, meshFaces);

            // ── Set outputs ───────────────────────────────────────────────────────────────
            dataAccess.SetListData(0, majorStressVectors);
            dataAccess.SetListData(1, minorStressVectors);
            dataAccess.SetListData(2, nodeVonMisesStress.Values.ToList());
            dataAccess.SetListData(3, originalVertices);
            dataAccess.SetListData(4, deformedVertices);
            dataAccess.SetData(5, deformedMesh);
        }
    }
}
