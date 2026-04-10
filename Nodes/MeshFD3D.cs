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
    // Unique identifier for this node. Never change this after the node is in use —
    // Synera uses it to reconnect saved graphs to the correct node class.
    [Guid("9D4E6F7A-0B1C-8D9E-0F1A-2B3C4D5E6F7A")]
    public sealed class MeshFD3D : Node
    {
        // The constructor defines the node name, description, and all input/output ports.
        public MeshFD3D() : base("Mesh FD 3D")
        {
            Category = WobisCategories.Aida;
            Subcategory = WobisSubcategories.Field;
            Description = "Deforms a shell mesh (quad or triangular) guided by both von Mises stress " +
                          "magnitude and principal stress directions, while keeping all nodes projected " +
                          "onto a 3D reference surface. Without surface projection, stress-guided " +
                          "relaxation on curved 3D geometries causes nodes to drift off the surface " +
                          "as the weighted averages are computed in free 3D space. This node corrects " +
                          "that by snapping each moved node back to the nearest point on the surface " +
                          "after every displacement step.";
            GuiPriority = 2;
            Keywords = "deform mesh;stress-guided;principal stress;anisotropic;3D surface;projection";

            // ── Input ports ───────────────────────────────────────────────────────────────

            // Port index 0: the solved FEA model providing stress results and mesh topology.
            InputParameterManager.AddParameter<IModel>("Input model",
                "Solved FEA model with shell elements whose stress results will guide the deformation.",
                ParameterAccess.Item);

            // Port index 1: the 3D reference surface onto which all nodes are projected after each move.
            // This is the key additional input compared to MeshFD.
            // The surface should closely match the geometry of the FEA mesh — ideally it is the
            // exact CAD surface from which the mesh was originally generated.
            InputParameterManager.AddParameter<IMesh>("Target surface",
                "3D reference surface mesh. After each relaxation step, every interior node is snapped " +
                "to the closest point on this surface to prevent drift off the geometry.",
                ParameterAccess.Item);

            // Port index 2: number of relaxation passes.
            InputParameterManager.AddParameter<SyneraInt>("Iterations",
                "Number of relaxation iterations. Higher values produce smoother, more stress-aligned meshes.",
                ParameterAccess.Item);

            // Port index 3: how strongly von Mises stress pulls nodes toward stressed regions.
            // Alpha = 0 → equal weighting (standard Laplacian smoothing).
            InputParameterManager.AddParameter<SyneraInt>("Alpha",
                "Controls how strongly von Mises stress influences node movement. " +
                "Alpha = 0 gives uniform Laplacian smoothing; higher values pull nodes toward high-stress regions.",
                ParameterAccess.Item);

            // Port index 4: how strongly principal stress directions steer node movement.
            // Beta = 0 → isotropic smoothing. Higher Beta → elements elongate along stress trajectories.
            // Step size is always preserved — only direction changes — so inversion is not possible.
            InputParameterManager.AddParameter<SyneraDouble>("Beta",
                "Controls how strongly principal stress directions influence node movement. " +
                "Beta = 0 disables directional anisotropy; higher values rotate the displacement " +
                "toward the principal stress axes without changing the step size.",
                ParameterAccess.Item);

            // Port index 5: damping factor for triangular meshes.
            // Limits how far each node moves per iteration to prevent sliver elements.
            // Lambda = 1.0 → full step (no damping). Lambda = 0.5 → half step.
            // Ignored for quad meshes (always treated as 1.0 internally).
            // Defaults to 1.0 so the port can be left unconnected for quad meshes.
            InputParameterManager.AddParameter<SyneraDouble>("Lambda",
                "Damping factor for triangular meshes — controls what fraction of the computed " +
                "displacement is applied each iteration. Valid range: (0, 1]. " +
                "Lambda = 1.0 applies the full step; Lambda = 0.5 applies half the step. " +
                "Ignored for quad meshes (always treated as 1.0). Defaults to 1.0.",
                ParameterAccess.Item,
                new SyneraDouble(1.0));

            // ── Output ports ──────────────────────────────────────────────────────────────

            // Port index 0: major principal stress vector per element.
            OutputParameterManager.AddParameter<Vector3D>("Max Principle Stress",
                "Maximum principal stress vectors per element.",
                ParameterAccess.List);

            // Port index 1: minor principal stress vector per element.
            OutputParameterManager.AddParameter<Vector3D>("Min Principle Stress",
                "Minimum principal stress vectors per element.",
                ParameterAccess.List);

            // Port index 2: von Mises stress averaged to each node from its surrounding elements.
            OutputParameterManager.AddParameter<SyneraDouble>("Von-Mises Stress",
                "Von Mises stress value averaged per node.",
                ParameterAccess.List);

            // Port index 3: original node positions before any deformation.
            OutputParameterManager.AddParameter<Point3D>("Vertices",
                "Original node locations before deformation.",
                ParameterAccess.List);

            // Port index 4: node positions after stress-guided relaxation and surface projection.
            OutputParameterManager.AddParameter<Point3D>("Deformed Vertices",
                "Node locations after stress-guided deformation, projected onto the target surface.",
                ParameterAccess.List);

            // Port index 5: the rebuilt mesh using the deformed, surface-projected vertex positions.
            OutputParameterManager.AddParameter<IMesh>("Mesh",
                "Resulting deformed mesh aligned with the principal stress field and lying on the target surface.",
                ParameterAccess.Item);
        }

        /// <summary>
        /// Converts a 6-component stress tensor into a single von Mises stress scalar.
        ///
        /// Von Mises stress combines normal and shear components into one value representing
        /// the overall stress intensity at a point. Used to identify high-stress regions
        /// where the mesh should be denser (pulled toward by the Alpha weighting).
        ///
        /// Input order: [Sxx, Syy, Szz, Sxy, Syz, Sxz]
        /// </summary>
        private static double CalculateVonMisesStress(double[] stressTensor)
        {
            double sxx = stressTensor[0], syy = stressTensor[1], szz = stressTensor[2];
            double sxy = stressTensor[3], syz = stressTensor[4], sxz = stressTensor[5];

            // Normal stress differences — how unevenly the material is stretched along each axis.
            double normalPart = (sxx - syy) * (sxx - syy)
                              + (syy - szz) * (syy - szz)
                              + (szz - sxx) * (szz - sxx);

            // Shear stress — how much the material is being twisted.
            double shearPart = 6.0 * (sxy * sxy + syz * syz + sxz * sxz);

            return Math.Sqrt(0.5 * (normalPart + shearPart));
        }

        /// <summary>
        /// Moves interior mesh nodes iteratively to align with the stress field,
        /// then projects each moved node back onto the 3D target surface.
        ///
        /// WHY SURFACE PROJECTION IS NEEDED:
        ///   On a flat (2D) mesh, all nodes already lie in the same plane, so weighted
        ///   averages of neighbor positions naturally stay in that plane. On a curved 3D
        ///   surface, this no longer holds: the weighted average of points on a curved
        ///   surface is a point INSIDE the surface (below it), not ON it. Without
        ///   correction, nodes gradually sink below the geometry with each iteration.
        ///
        ///   The fix is simple: after computing the new position for a node, find the
        ///   closest point on the reference surface and snap the node there. This "project
        ///   and correct" step runs after every displacement, keeping all interior nodes
        ///   on the surface throughout the relaxation.
        ///
        /// UPDATE SCHEME — controlled by useJacobi:
        ///   Gauss-Seidel (useJacobi = false) — quad / mixed meshes:
        ///     Positions updated in-place immediately. Subsequent nodes see updated neighbors.
        ///     Fast convergence. Stable for ~4 neighbors at 90°.
        ///
        ///   Jacobi (useJacobi = true) — triangular meshes:
        ///     All new positions computed from a frozen snapshot, applied simultaneously.
        ///     Prevents chaotic propagation through dense triangular connectivity (~6 neighbors at 60°).
        ///
        /// ALGORITHM PER NODE PER ITERATION:
        ///   1. Compute a von Mises weighted centroid of neighbor positions.
        ///   2. Decompose the displacement (node → centroid) into major/minor stress components.
        ///   3. Amplify stress-aligned components by Beta, rescale to original length (no overshoot).
        ///   4. Apply Lambda fraction of the displacement (damping for triangular meshes).
        ///   5. Project the candidate position onto the target surface. ← unique to MeshFD3D
        ///   6. Store or apply the projected position.
        /// </summary>
        private static void DeformNodes(
            List<INode> allNodes,           // All nodes in the FEA mesh
            int iterations,                 // Number of relaxation passes
            double alpha,                   // Von Mises stress influence strength
            double beta,                    // Principal stress direction bias strength
            double lambda,                  // Damping factor — fraction of displacement applied per iteration
            bool useJacobi,                 // true = Jacobi (tri), false = Gauss-Seidel (quad/mixed)
            IMesh targetSurface,            // The 3D surface nodes are projected onto after each move
            Dictionary<int, double> nodeVonMisesStress,
            Dictionary<int, List<int>> nodeToNeighborsMap,
            HashSet<int> boundaryNodeIndices,
            Dictionary<int, (Vector3D majorDirection, Vector3D minorDirection)> nodePrincipalDirections)
        {
            // Normalise von Mises stress to [0, 1] so Alpha is unit-independent.
            double minStress = nodeVonMisesStress.Values.Min();
            double maxStress = nodeVonMisesStress.Values.Max();
            double stressRange = maxStress - minStress + 1e-6; // epsilon prevents division by zero

            // nodeDict: live node lookup for Gauss-Seidel (reflects in-place updates immediately).
            Dictionary<int, INode> nodeDict = allNodes.ToDictionary(n => n.Index);

            // currentPositions: frozen snapshot for Jacobi (never updated mid-iteration).
            // All centroid computations within one Jacobi iteration read from this dictionary.
            Dictionary<int, Point3D> currentPositions = allNodes.ToDictionary(n => n.Index, n => n.Location);

            for (int iteration = 0; iteration < iterations; iteration++)
            {
                // Jacobi write buffer — only used when useJacobi is true.
                // Stores the new projected positions for all nodes before applying them.
                Dictionary<int, Point3D> newPositions = useJacobi
                    ? new Dictionary<int, Point3D>(allNodes.Count)
                    : null;

                foreach (INode node in allNodes)
                {
                    // Boundary nodes stay fixed — they define the mesh outline and are already
                    // on the surface by definition. No projection needed for them.
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

                    // Retrieve the averaged principal stress directions for this node.
                    // Computed earlier by averaging element-level stress over all connected elements.
                    nodePrincipalDirections.TryGetValue(node.Index, out var principalDirs);
                    Vector3D majorDir = principalDirs.majorDirection;
                    Vector3D minorDir = principalDirs.minorDirection;

                    double majorMagnitude = majorDir.Length;
                    double minorMagnitude = minorDir.Length;

                    // Per-node local normalisation: the dominant stress direction always scores 1.0.
                    // This ensures Beta has consistent influence regardless of absolute stress magnitudes.
                    double localMax = Math.Max(majorMagnitude, minorMagnitude) + 1e-6;
                    double normalizedMajorMagnitude = majorMagnitude / localMax;
                    double normalizedMinorMagnitude = minorMagnitude / localMax;

                    // Unit vectors along each principal stress direction.
                    // Fall back to world axes if the node has negligible stress.
                    Vector3D majorUnit = majorMagnitude > 1e-10 ? majorDir.Normalized() : new Vector3D(1, 0, 0);
                    Vector3D minorUnit = minorMagnitude > 1e-10 ? minorDir.Normalized() : new Vector3D(0, 1, 0);

                    // ── Step 1: Von Mises weighted Laplacian centroid ─────────────────────────
                    // Compute a weighted average of neighbor positions.
                    // Gauss-Seidel reads live positions from nodeDict (some already updated this pass).
                    // Jacobi reads from the frozen snapshot in currentPositions (all from start of pass).
                    Point3D weightedPositionSum = new Point3D(0, 0, 0);
                    double totalWeight = 0.0;

                    foreach (int neighborIdx in neighbors)
                    {
                        Point3D neighborPos = useJacobi
                            ? currentPositions[neighborIdx]
                            : nodeDict[neighborIdx].Location;

                        // Neighbors with higher von Mises stress get a larger weight.
                        // At Alpha = 0, all neighbors have equal weight = 1.0 (pure Laplacian).
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

                    // This is where the node would land under pure von Mises weighted smoothing.
                    // NOTE: on a curved 3D surface, this centroid already drifts slightly off the
                    // surface. The surface projection in Step 5 corrects this.
                    Point3D centroid = weightedPositionSum / totalWeight;

                    // ── Step 2: Decompose displacement into principal stress components ─────────
                    Point3D oldPosition = useJacobi ? currentPositions[node.Index] : node.Location;
                    Vector3D delta = Point3D.Subtract(centroid, oldPosition);

                    // Project delta onto each principal stress direction.
                    double dMajor = delta * majorUnit;
                    double dMinor = delta * minorUnit;

                    // The residual is the part of delta not aligned with either stress direction.
                    // On a curved surface this includes curvature-induced out-of-plane motion,
                    // which the surface projection step will correct.
                    Vector3D deltaRest = delta - dMajor * majorUnit - dMinor * minorUnit;

                    // ── Step 3: Steer direction toward principal stress axes ───────────────────
                    // Beta amplifies the stress-aligned components. The vector is then rescaled
                    // back to the original delta length so Beta changes direction only, not magnitude.
                    // This prevents overshoot and mesh inversion regardless of how large Beta is.
                    Vector3D anisotropicDirection = deltaRest
                        + (1.0 + beta * normalizedMajorMagnitude) * dMajor * majorUnit
                        + (1.0 + beta * normalizedMinorMagnitude) * dMinor * minorUnit;

                    double deltaLength = delta.Length;
                    double anisotropicLength = anisotropicDirection.Length;

                    Vector3D anisotropicDelta;
                    if (anisotropicLength > 1e-10 && deltaLength > 1e-10)
                        anisotropicDelta = (deltaLength / anisotropicLength) * anisotropicDirection;
                    else
                        anisotropicDelta = delta; // fallback: use the unsteered displacement

                    // ── Step 4: Apply damping ─────────────────────────────────────────────────
                    // Lambda < 1.0 limits per-iteration displacement to prevent sliver elements
                    // in triangular meshes. For quads, lambda is always 1.0.
                    Point3D candidatePosition = oldPosition + lambda * anisotropicDelta;

                    // ── Step 5: Project candidate position onto the target surface ─────────────
                    // This is the step that makes MeshFD3D different from MeshFD.
                    //
                    // WHY: In 2D (flat) meshes, weighted averages of surface points stay on the
                    // surface naturally. On a curved 3D surface, the weighted average of points
                    // on the surface lies BELOW the surface (inside the curvature). Each iteration
                    // pushes nodes slightly inward. After many iterations the mesh visibly sinks.
                    //
                    // HOW: Find the closest point on the reference surface mesh to the candidate
                    // position. Move the node there instead of to the raw candidate position.
                    // This "lift" step restores the node onto the surface after each displacement.
                    //
                    // BOUNDARY NODES: they are skipped entirely (see top of loop), so they are
                    // never moved and never need projection. They remain exactly where they started.
                    double distanceToSurface;
                    Point3D projectedPosition = targetSurface.ClosestPoint(candidatePosition, out distanceToSurface);

                    // ── Step 6: Store or apply the projected position ─────────────────────────
                    if (useJacobi)
                        newPositions[node.Index] = projectedPosition;
                    else
                        node.Location = projectedPosition;
                }

                // ── Jacobi write: apply all projected positions simultaneously ──────────────────
                // Only after ALL nodes have been computed and projected do we update positions.
                // No node in this iteration influenced another node's result.
                if (useJacobi)
                {
                    foreach (INode node in allNodes)
                    {
                        if (newPositions.TryGetValue(node.Index, out Point3D newPos))
                        {
                            node.Location = newPos;
                            // Update the read buffer so the next iteration starts from
                            // the positions we just applied.
                            currentPositions[node.Index] = newPos;
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Called by Synera whenever the node needs to recompute.
        /// Reads and validates inputs, detects mesh type, runs relaxation with
        /// surface projection, and writes outputs.
        /// </summary>
        protected override void SolveInstance(IDataAccess dataAccess)
        {
            // ── Read all inputs from connected ports ──────────────────────────────────────
            bool isDataSuccess = dataAccess.GetData(0, out IModel model);
            isDataSuccess &= dataAccess.GetData(1, out IMesh targetSurface);
            isDataSuccess &= dataAccess.GetData(2, out SyneraInt iterations);
            isDataSuccess &= dataAccess.GetData(3, out SyneraInt alpha);
            isDataSuccess &= dataAccess.GetData(4, out double beta);
            isDataSuccess &= dataAccess.GetData(5, out double lambda);

            // If any required port is disconnected or has no data, stop here.
            if (!isDataSuccess)
                return;

            // ── Validate all inputs ───────────────────────────────────────────────────────
            // Collect all errors before returning so the user sees everything at once.
            bool isValid = true;

            if (!model.IsValid)
            {
                AddError(0, "Invalid model.");
                isValid = false;
            }
            else if (!model.HasShellElementsOnly())
            {
                // Shell elements define the surface mesh we are deforming.
                // Solid elements have 3D topology that is not compatible with this approach.
                AddError(0, "A model containing shell elements only is expected.");
                isValid = false;
            }

            if (targetSurface == null)
            {
                AddError(1, "A target surface mesh must be connected.");
                isValid = false;
            }

            if (iterations <= 0)
            {
                AddError(2, "Iterations must be a positive integer.");
                isValid = false;
            }

            if (alpha < 0)
            {
                AddError(3, "Alpha must be a non-negative integer.");
                isValid = false;
            }

            if (double.IsPositiveInfinity(alpha))
            {
                AddError(3, "Alpha is infinite. Please provide a finite value.");
                isValid = false;
            }

            if (beta < 0)
            {
                AddError(4, "Beta must be a non-negative value.");
                isValid = false;
            }

            if (double.IsPositiveInfinity(beta))
            {
                AddError(4, "Beta is infinite. Please provide a finite value.");
                isValid = false;
            }

            if (lambda <= 0 || lambda > 1.0)
            {
                AddError(5, "Lambda must be in the range (0, 1]. Use 0.5 as a safe starting value for triangular meshes.");
                isValid = false;
            }

            if (!isValid)
                return;

            // ── Detect mesh type and configure update strategy ────────────────────────────
            // Pure triangular → Jacobi update + Lambda damping.
            // Quad or mixed   → Gauss-Seidel update, Lambda forced to 1.0.
            // Mixed meshes are treated as quad-dominant for stability.
            IEnumerable<IShellElement> shellElements = model.Elements.OfType<IShellElement>();

            bool isTriangular = shellElements.All(elem =>
                elem.NodeIndices.Except(elem.GetMidNodeIndices()).Count() == 3);

            bool isMixed = !isTriangular && shellElements.Any(elem =>
                elem.NodeIndices.Except(elem.GetMidNodeIndices()).Count() == 3);

            bool useJacobi = isTriangular;
            double effectiveLambda = isTriangular ? lambda : 1.0;

            if (isTriangular)
                AddNotification("Triangular mesh detected — Jacobi update with Lambda damping applied.");
            else if (isMixed)
                AddNotification("Mixed quad/triangular mesh detected — treated as quad-dominant. Gauss-Seidel update applied. Lambda ignored.");
            else
                AddNotification("Quad mesh detected — Gauss-Seidel update applied. Lambda ignored.");

            // ── Step 1: Compute principal stress vectors per element ──────────────────────
            // StressStreamlineBuilder performs eigenvalue decomposition on each element's
            // stress tensor to extract the two in-plane principal stress directions.
            Progress progress = this.CreateProgress();
            StressStreamlineBuilder streamlineBuilder = new StressStreamlineBuilder(model);

            Progress stressProgress = progress.CreateSubtask(0.5);
            Dictionary<int, (Vector3D majorVector, Vector3D minorVector)> elementPrincipalStress =
                streamlineBuilder.GetPrincipalStressValues(stressProgress);

            List<Vector3D> majorStressVectors = elementPrincipalStress.Values.Select(v => v.majorVector).ToList();
            List<Vector3D> minorStressVectors  = elementPrincipalStress.Values.Select(v => v.minorVector).ToList();

            // ── Step 2: Build mesh topology lookup maps ───────────────────────────────────
            List<INode> nodes = model.Nodes.ToList();

            // Map: node index → indices of all elements sharing that node.
            // Only corner nodes are included (mid-side nodes of higher-order elements excluded).
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

            // Map: node index → indices of all directly connected (edge-adjacent) neighbor nodes.
            // Defines the "neighborhood" used in the Laplacian relaxation.
            Dictionary<int, List<int>> nodeToNeighborsMap = new Dictionary<int, List<int>>();
            foreach (IndexPair edge in model.Edges)
            {
                int a = edge.A;
                int b = edge.B;

                if (a == b) // skip degenerate zero-length edges
                    continue;

                if (!nodeToNeighborsMap.ContainsKey(a)) nodeToNeighborsMap[a] = new List<int>();
                if (!nodeToNeighborsMap.ContainsKey(b)) nodeToNeighborsMap[b] = new List<int>();

                // Each edge is undirected — add each node as the other's neighbor.
                if (!nodeToNeighborsMap[a].Contains(b)) nodeToNeighborsMap[a].Add(b);
                if (!nodeToNeighborsMap[b].Contains(a)) nodeToNeighborsMap[b].Add(a);
            }

            // ── Step 3: Compute per-node von Mises stress ─────────────────────────────────
            // FEA provides stress per element. Average it over connected elements to get
            // a smooth scalar stress field defined at each node.
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
                    // Nodes with no connected elements (isolated) carry zero stress.
                    nodeVonMisesStress[node.Index] = 0.0;
                }
            }

            // ── Step 4: Average principal stress directions per node ───────────────────────
            // Interpolate element-level principal stress vectors to the nodes by averaging
            // over all elements connected to each node.
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

            // ── Step 5: Identify fixed boundary nodes ─────────────────────────────────────
            // Nodes on naked (free/open boundary) edges define the mesh outline.
            // They are kept fixed to preserve the footprint of the structure.
            // Because they never move, they never drift off the surface — no projection needed.
            IEnumerable<int> nakedEdgeIndices = streamlineBuilder.Mesh.GetNakedEdgeIndices();
            HashSet<int> boundaryNodeIndices = new HashSet<int>();
            foreach (int edgeIdx in nakedEdgeIndices)
            {
                IndexPair edge = model.GetNodeIndicesForEdge(edgeIdx);
                boundaryNodeIndices.AddRange(edge.ToList());
            }

            // Record original positions before any deformation.
            List<Point3D> originalVertices = nodes.Select(n => n.Location).ToList();

            // ── Step 6: Run stress-guided relaxation with surface projection ───────────────
            // This is where the mesh deformation actually happens.
            // Each interior node is moved toward a stress-weighted centroid of its neighbors,
            // steered toward the principal stress axes by Beta, damped by Lambda,
            // and then snapped back onto the target surface.
            DeformNodes(
                nodes,
                iterations,
                alpha,
                beta,
                effectiveLambda,
                useJacobi,
                targetSurface,
                nodeVonMisesStress,
                nodeToNeighborsMap,
                boundaryNodeIndices,
                nodePrincipalDirections);

            List<Point3D> deformedVertices = nodes.Select(n => n.Location).ToList();

            // ── Step 7: Rebuild the mesh from the deformed, projected vertex positions ──────
            // Mesh connectivity (face definitions) is unchanged — only positions differ.
            List<MeshFace> meshFaces = shellElements
                .Select(elem => new MeshFace(elem.NodeIndices.Select(i => i - 1).ToArray()))
                .ToList();

            IMeshKernel meshKernel = Application.Current.KernelManager.Get<IMeshKernel>();
            IMesh deformedMesh = meshKernel.CreateFromVerticesAndFaces(deformedVertices, meshFaces);

            // ── Step 8: Write all outputs ─────────────────────────────────────────────────
            dataAccess.SetListData(0, majorStressVectors);
            dataAccess.SetListData(1, minorStressVectors);
            dataAccess.SetListData(2, nodeVonMisesStress.Values.ToList());
            dataAccess.SetListData(3, originalVertices);
            dataAccess.SetListData(4, deformedVertices);
            dataAccess.SetData(5, deformedMesh);
        }
    }
}
