using System;
using System.Linq;
using System.Collections.Generic;
using Synera.Core.Implementation.ApplicationService;
using Synera.Kernels;
using Synera.Kernels.DataTypes;
using Synera.Kernels.Fem.Elements;
using Synera.Kernels.Fem.Model;
using Synera.Kernels.Fem.Results;
using Synera.Kernels.Mesh;
using Synera.Utilities;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using MathNet.Numerics.LinearAlgebra.Factorization;

namespace AidaTool.DataTypes
{
    public class StressStreamlineBuilder
    {
        public IModel Model { get; }

        public IMesh Mesh { get; private set; }

        public List<int> FaceToElementIndices { get; private set; }

        public HashSet<int> OpenMeshEdges { get; private set; }

        public StressStreamlineBuilder(IModel model)
        {
            Model = model;
            List<int> meshFaceToElementIndexMap;
            Mesh = ConvertShellsToMesh(model, out meshFaceToElementIndexMap);
            FaceToElementIndices = meshFaceToElementIndexMap;
            OpenMeshEdges = new HashSet<int>(Mesh.GetNakedEdgeIndices());
        }

        private static IMesh ConvertShellsToMesh(IModel model, out List<int> meshFaceToElementIndexMap)
        {
            // Collection to store the unique 3D point coordinates that will form the vertices of the new mesh.
            List<Point3D> meshVertices = new List<Point3D>(model.NodeCount);

            // Collection to store the definitions of each face (polygon) in the new mesh.
            List<MeshFace> meshFaces = new List<MeshFace>(model.ElementCount);

            // A dictionary to efficiently map an original FEA node ID to its new, compact vertex index
            // within the 'meshVertices' list. This prevents duplicate vertices.
            Dictionary<int, int> originalNodeIdToNewVertexIndexMap = new Dictionary<int, int>(model.NodeCount);

            // Initialize the output list that will store the mapping from new mesh face indices
            // back to their corresponding original FEA element IDs.
            meshFaceToElementIndexMap = new List<int>(model.ElementCount);

            // Iterate through each IShellElement present in the FEA model.
            foreach (IShellElement shellElement in model.Elements.OfType<IShellElement>())
            {
                // List to temporarily hold the new vertex indices for the current mesh face being constructed.
                // Shell elements are typically 3-noded (tri) or 4-noded (quad).
                List<int> currentFaceVertexIndices = new List<int>(4);

                // Extract only the corner node IDs from the current shell element.
                // Higher-order elements (e.g., 8-noded quads) will be simplified to their 4-noded base.
                IEnumerable<int> cornerNodeIds = shellElement.NodeIndices.Except(shellElement.GetMidNodeIndices());

                foreach (int originalNodeId in cornerNodeIds)
                {
                    int newVertexIndex;

                    // Check if this original FEA node has already been added to our mesh's vertex list.
                    if (!originalNodeIdToNewVertexIndexMap.TryGetValue(originalNodeId, out newVertexIndex))
                    {
                        // If it's a new node (not yet mapped):
                        // 1. Assign it a new, sequential index in our 'meshVertices' list.
                        newVertexIndex = meshVertices.Count;
                        // 2. Store this mapping: original FEA node ID -> new vertex index.
                        originalNodeIdToNewVertexIndexMap.Add(originalNodeId, newVertexIndex);
                        // 3. Add the node's geometric location to our list of mesh vertices.
                        meshVertices.Add(model.GetNode(originalNodeId).Location);
                    }
                    // Add the determined vertex index (either new or existing) to the current face's vertex list.
                    currentFaceVertexIndices.Add(newVertexIndex);
                }

                // After gathering all corner vertices for the current shell element:
                // 1. Add the original FEA element's ID to the mapping list. The position in this list
                //    corresponds to the index of the MeshFace that will be created next.
                meshFaceToElementIndexMap.Add(shellElement.Index);
                // 2. Create a new MeshFace object from the collected vertex indices and add it
                //    to the list that defines the faces of our new geometric mesh.
                meshFaces.Add(new MeshFace(currentFaceVertexIndices.ToArray()));
            }

            // Finally, use the Synera MeshKernel to construct and return the complete IMesh
            // from the assembled unique vertices and face definitions.
            var meshKernel = Application.Current.KernelManager.Get<IMeshKernel>();
            var mesh = meshKernel.CreateFromVerticesAndFaces(meshVertices, meshFaces);
            return mesh;
        }

        public Dictionary<int, (Vector3D major, Vector3D minor)> GetPrincipalStressValues(Progress progress)
        {
            // Retrieve the raw element stress tensor results from the model.
            // This assumes the stress tensor results are available and correctly named.
            IResultsAtTimeStep elementStressResults = Model.Results
                .First(result => result.Name.Equals(ResultType.ElementStressTensor.Name))
                .Values.First() // Get the first result set (e.g., for a single load case)
                .Values.First(); // Get the results for the first time step

            // Dictionary to store the final calculated principal stress directions for each element.
            // Key: Original FEA Element ID
            // Value: Tuple containing (Major Principal Stress Vector, Minor Principal Stress Vector)
            Dictionary<int, (Vector3D majorDirection, Vector3D minorDirection)> principalStressDirections =
                new Dictionary<int, (Vector3D, Vector3D)>(elementStressResults.Count);

            // Counter for tracking processed elements to report progress.
            int elementsProcessedCount = 0;

            // Iterate through each element's stress tensor data.
            // Each entry contains the element ID and an array of stress components.
            foreach (KeyValuePair<int, double[]> elementStressData in elementStressResults)
            {
                // The raw stress data might contain values from multiple integration points.
                // We need to average these components to get a single, representative stress tensor for the element.
                // The stress tensor is typically represented by 6 unique components: [Sxx, Syy, Szz, Sxy, Syz, Sxz].
                double[] averagedStressComponents = new double[6];
                int componentsPerPoint = 6; // Number of stress components (Sxx, Syy, Szz, Sxy, Syz, Sxz)
                int numberOfIntegrationPoints = elementStressData.Value.Length / (componentsPerPoint * 3); // Assuming 3 axes per component for full tensor

                for (int i = 0; i < componentsPerPoint; i++)
                {
                    // Average each stress component across all integration points.
                    averagedStressComponents[i] = Enumerable.Range(0, numberOfIntegrationPoints)
                        .Select(pointIdx => elementStressData.Value[componentsPerPoint * pointIdx + i])
                        .Average();
                }

                // Construct the 3x3 stress tensor matrix from the averaged components.
                // MathNet.Numerics.LinearAlgebra expects column-major order by default.
                // The matrix is symmetric:
                // [ Sxx  Sxy  Sxz ]
                // [ Sxy  Syy  Syz ]
                // [ Sxz  Syz  Szz ]
                Matrix<double> stressTensorMatrix = new DenseMatrix(3, 3);
                stressTensorMatrix[0, 0] = averagedStressComponents[0]; // Sxx
                stressTensorMatrix[1, 0] = averagedStressComponents[3]; // Sxy
                stressTensorMatrix[2, 0] = averagedStressComponents[5]; // Sxz

                stressTensorMatrix[0, 1] = averagedStressComponents[3]; // Syx (same as Sxy)
                stressTensorMatrix[1, 1] = averagedStressComponents[1]; // Syy
                stressTensorMatrix[2, 1] = averagedStressComponents[4]; // Syz

                stressTensorMatrix[0, 2] = averagedStressComponents[5]; // Szx (same as Sxz)
                stressTensorMatrix[1, 2] = averagedStressComponents[4]; // Szy (same as Syz)
                stressTensorMatrix[2, 2] = averagedStressComponents[2]; // Szz

                // Perform Eigenvalue Decomposition (EVD) on the stress tensor matrix.
                // For a symmetric matrix, EVD yields real eigenvalues (principal stresses)
                // and orthogonal eigenvectors (principal stress directions).
                Evd<double> eigenvalueDecomposition = stressTensorMatrix.Evd(Symmetricity.Symmetric);

                // For shell elements, we are interested in the in-plane principal stress directions.
                // One of the three eigenvectors will typically be perpendicular to the shell's surface.
                // We identify and exclude this "through-thickness" eigenvector.

                // Get the normal vector of the mesh face corresponding to the current element.
                // Assuming a 1-to-1 mapping where mesh face index is (element ID - 1).
                Vector3D meshFaceNormal = Mesh.GetFaceNormal(elementStressData.Key - 1);

                // Create a list of indices for the three eigenvectors (0, 1, 2).
                List<int> eigenvectorIndexList = Enumerable.Range(0, 3).ToList();

                // Calculate the absolute dot product of each eigenvector with the mesh face normal.
                // The eigenvector with the largest absolute dot product is most aligned with the normal.
                double[] dotProductsWithNormal = new double[3];
                for (int i = 0; i < 3; i++)
                {
                    Vector3D currentEigenvector = new Vector3D(
                        eigenvalueDecomposition.EigenVectors.At(0, i),
                        eigenvalueDecomposition.EigenVectors.At(1, i),
                        eigenvalueDecomposition.EigenVectors.At(2, i));
                    dotProductsWithNormal[i] = Math.Abs(currentEigenvector * meshFaceNormal);
                }

                // Find the index of the eigenvector that is most perpendicular to the shell plane (aligned with normal).
                int throughThicknessVectorIndex = Array.IndexOf(dotProductsWithNormal, dotProductsWithNormal.Max());

                // Remove this eigenvector's index from our list. The remaining two indices correspond
                // to the in-plane principal stress directions.
                eigenvectorIndexList.Remove(throughThicknessVectorIndex);

                // Extract and normalize the two in-plane principal stress vectors.
                // The order (major/minor) is typically determined by the magnitude of their corresponding eigenvalues,
                // but here we simply take the two remaining directions.
                Vector3D firstInPlaneVector = new Vector3D(
                    eigenvalueDecomposition.EigenVectors.At(0, eigenvectorIndexList[0]),
                    eigenvalueDecomposition.EigenVectors.At(1, eigenvectorIndexList[0]),
                    eigenvalueDecomposition.EigenVectors.At(2, eigenvectorIndexList[0])).Normalized();

                Vector3D secondInPlaneVector = new Vector3D(
                    eigenvalueDecomposition.EigenVectors.At(0, eigenvectorIndexList[1]),
                    eigenvalueDecomposition.EigenVectors.At(1, eigenvectorIndexList[1]),
                    eigenvalueDecomposition.EigenVectors.At(2, eigenvectorIndexList[1])).Normalized();

                // Store the calculated principal stress directions for the current element.
                // The assignment (first/second to major/minor) depends on the eigenvalue order,
                // but for streamline tracing, the two orthogonal in-plane directions are what's needed.
                //principalStressDirections[elementStressData.Key] = (firstInPlaneVector, secondInPlaneVector);

                // NEW CODE 
                // Eigen values
                Vector<double> eigenValues = eigenvalueDecomposition.EigenValues.Real();
                principalStressDirections[elementStressData.Key] = (eigenValues.At(eigenvectorIndexList[0]) * firstInPlaneVector, eigenValues.At(eigenvectorIndexList[1]) * secondInPlaneVector);
                // NEW CODE

                elementsProcessedCount++;
                // Report progress to the caller, updating every 4092 elements (arbitrary batch size).
                if (elementsProcessedCount % 4092 == 0)
                {
                    progress.SetProgress(elementsProcessedCount / (double)elementStressResults.Count);
                }
            }

            // Return the dictionary containing all calculated principal stress directions.
            return principalStressDirections;
        }

        //public Dictionary<int, (Vector3D major, Vector3D minor)> GetPrincipalStressVectors(Progress progress)
        //{
        //    // Retrieve the raw element stress tensor results from the model.
        //    // This assumes the stress tensor results are available and correctly named.
        //    IResultsAtTimeStep elementStressResults = Model.Results
        //        .First(result => result.Name.Equals(ResultType.ElementStressTensor.Name))
        //        .Values.First() // Get the first result set (e.g., for a single load case)
        //        .Values.First(); // Get the results for the first time step

        //    // Dictionary to store the final calculated principal stress directions for each element.
        //    // Key: Original FEA Element ID
        //    // Value: Tuple containing (Major Principal Stress Vector, Minor Principal Stress Vector)
        //    Dictionary<int, (Vector3D majorDirection, Vector3D minorDirection)> principalStressDirections =
        //        new Dictionary<int, (Vector3D, Vector3D)>(elementStressResults.Count);

        //    // Counter for tracking processed elements to report progress.
        //    int elementsProcessedCount = 0;

        //    // Iterate through each element's stress tensor data.
        //    // Each entry contains the element ID and an array of stress components.
        //    foreach (KeyValuePair<int, double[]> elementStressData in elementStressResults)
        //    {
        //        // The raw stress data might contain values from multiple integration points.
        //        // We need to average these components to get a single, representative stress tensor for the element.
        //        // The stress tensor is typically represented by 6 unique components: [Sxx, Syy, Szz, Sxy, Syz, Sxz].
        //        double[] averagedStressComponents = new double[6];
        //        int componentsPerPoint = 6; // Number of stress components (Sxx, Syy, Szz, Sxy, Syz, Sxz)
        //        int numberOfIntegrationPoints = elementStressData.Value.Length / (componentsPerPoint * 3); // Assuming 3 axes per component for full tensor

        //        for (int i = 0; i < componentsPerPoint; i++)
        //        {
        //            // Average each stress component across all integration points.
        //            averagedStressComponents[i] = Enumerable.Range(0, numberOfIntegrationPoints)
        //                .Select(pointIdx => elementStressData.Value[componentsPerPoint * pointIdx + i])
        //                .Average();
        //        }

        //        // Construct the 3x3 stress tensor matrix from the averaged components.
        //        // MathNet.Numerics.LinearAlgebra expects column-major order by default.
        //        // The matrix is symmetric:
        //        // [ Sxx  Sxy  Sxz ]
        //        // [ Sxy  Syy  Syz ]
        //        // [ Sxz  Syz  Szz ]
        //        Matrix<double> stressTensorMatrix = new DenseMatrix(3, 3);
        //        stressTensorMatrix[0, 0] = averagedStressComponents[0]; // Sxx
        //        stressTensorMatrix[1, 0] = averagedStressComponents[3]; // Sxy
        //        stressTensorMatrix[2, 0] = averagedStressComponents[5]; // Sxz

        //        stressTensorMatrix[0, 1] = averagedStressComponents[3]; // Syx (same as Sxy)
        //        stressTensorMatrix[1, 1] = averagedStressComponents[1]; // Syy
        //        stressTensorMatrix[2, 1] = averagedStressComponents[4]; // Syz

        //        stressTensorMatrix[0, 2] = averagedStressComponents[5]; // Szx (same as Sxz)
        //        stressTensorMatrix[1, 2] = averagedStressComponents[4]; // Szy (same as Syz)
        //        stressTensorMatrix[2, 2] = averagedStressComponents[2]; // Szz

        //        // Perform Eigenvalue Decomposition (EVD) on the stress tensor matrix.
        //        // For a symmetric matrix, EVD yields real eigenvalues (principal stresses)
        //        // and orthogonal eigenvectors (principal stress directions).
        //        Evd<double> eigenvalueDecomposition = stressTensorMatrix.Evd(Symmetricity.Symmetric);

        //        // For shell elements, we are interested in the in-plane principal stress directions.
        //        // One of the three eigenvectors will typically be perpendicular to the shell's surface.
        //        // We identify and exclude this "through-thickness" eigenvector.

        //        // Get the normal vector of the mesh face corresponding to the current element.
        //        // Assuming a 1-to-1 mapping where mesh face index is (element ID - 1).
        //        Vector3D meshFaceNormal = Mesh.GetFaceNormal(elementStressData.Key - 1);

        //        // Create a list of indices for the three eigenvectors (0, 1, 2).
        //        List<int> eigenvectorIndexList = Enumerable.Range(0, 3).ToList();

        //        // Calculate the absolute dot product of each eigenvector with the mesh face normal.
        //        // The eigenvector with the largest absolute dot product is most aligned with the normal.
        //        double[] dotProductsWithNormal = new double[3];
        //        for (int i = 0; i < 3; i++)
        //        {
        //            Vector3D currentEigenvector = new Vector3D(
        //                eigenvalueDecomposition.EigenVectors.At(0, i),
        //                eigenvalueDecomposition.EigenVectors.At(1, i),
        //                eigenvalueDecomposition.EigenVectors.At(2, i));
        //            dotProductsWithNormal[i] = Math.Abs(currentEigenvector * meshFaceNormal);
        //        }

        //        // Find the index of the eigenvector that is most perpendicular to the shell plane (aligned with normal).
        //        int throughThicknessVectorIndex = Array.IndexOf(dotProductsWithNormal, dotProductsWithNormal.Max());

        //        // Remove this eigenvector's index from our list. The remaining two indices correspond
        //        // to the in-plane principal stress directions.
        //        eigenvectorIndexList.Remove(throughThicknessVectorIndex);

        //        // Extract and normalize the two in-plane principal stress vectors.
        //        // The order (major/minor) is typically determined by the magnitude of their corresponding eigenvalues,
        //        // but here we simply take the two remaining directions.
        //        Vector3D firstInPlaneVector = new Vector3D(
        //            eigenvalueDecomposition.EigenVectors.At(0, eigenvectorIndexList[0]),
        //            eigenvalueDecomposition.EigenVectors.At(1, eigenvectorIndexList[0]),
        //            eigenvalueDecomposition.EigenVectors.At(2, eigenvectorIndexList[0])).Normalized();

        //        Vector3D secondInPlaneVector = new Vector3D(
        //            eigenvalueDecomposition.EigenVectors.At(0, eigenvectorIndexList[1]),
        //            eigenvalueDecomposition.EigenVectors.At(1, eigenvectorIndexList[1]),
        //            eigenvalueDecomposition.EigenVectors.At(2, eigenvectorIndexList[1])).Normalized();

        //        // NEW CODE 
        //        // Eigen values
        //        var eigenValues = eigenvalueDecomposition.EigenValues.Real();

        //        // Store the calculated principal stress directions for the current element.
        //        // The assignment (first/second to major/minor) depends on the eigenvalue order,
        //        // but for streamline tracing, the two orthogonal in-plane directions are what's needed.
        //        //principalStressDirections[elementStressData.Key] = (firstInPlaneVector, secondInPlaneVector);
        //        principalStressDirections[elementStressData.Key] = (firstInPlaneVector, secondInPlaneVector);

        //        elementsProcessedCount++;
        //        // Report progress to the caller, updating every 4092 elements (arbitrary batch size).
        //        if (elementsProcessedCount % 4092 == 0)
        //        {
        //            progress.SetProgress(elementsProcessedCount / (double)elementStressResults.Count);
        //        }
        //    }

        //    // Return the dictionary containing all calculated principal stress directions.
        //    return principalStressDirections;
        //}

        //public ICurve CreateStreamline(IModel modelData, IDictionary<int, Vector3D> stressVectors, Point3D initialSeedPoint)
        //{
        //    double effectiveStepLength;
        //    // Calculate the average length of all edges in the FEA model.
        //    // Using half of this average as a step length often provides a good visual balance.
        //    if (modelData.Edges.Any())
        //    {
        //        effectiveStepLength = 0.5 * modelData.Edges.Average(edge =>
        //        {
        //            Point3D startNodeLocation = modelData.GetNode(edge.A).Location;
        //            Point3D endNodeLocation = modelData.GetNode(edge.B).Location;
        //            return startNodeLocation.DistanceTo(endNodeLocation);
        //        });
        //    }
        //    else
        //    {
        //        // Fallback to a small arbitrary value if the model has no edges (e.g., a point model).
        //        effectiveStepLength = 0.1;
        //    }

        //    // Calculate the maximum number of steps allowed for the streamline to prevent excessive length or infinite loops.
        //    // This is a heuristic based on the model's bounding box diagonal.
        //    BoxDef modelBoundingBox = Mesh.GetBoundingBox();
        //    double modelExtent = modelBoundingBox.Diagonal.Length;
        //    int maximumTraceSteps = 10 * (int)Math.Ceiling(modelExtent / effectiveStepLength);

        //    // Project the initial seed point onto the mesh surface to ensure the streamline starts directly on the geometry.
        //    double distanceToMesh; // This output parameter is not strictly needed for the logic here.
        //    Point3D startPointOnMesh = Mesh.ClosestPoint(initialSeedPoint, out distanceToMesh);

        //    // Trace the streamline segment in the positive direction of the stress vectors.
        //    List<Point3D> forwardPathSegment = TraceStreamlineSegment(
        //        startPointOnMesh,
        //        stressVectors,
        //        effectiveStepLength,
        //        maximumTraceSteps,
        //        true); // true for positive direction

        //    // Trace the streamline segment in the negative direction of the stress vectors.
        //    // Note: The TraceStreamlineSegment function is designed to handle this by using the start point
        //    // and then moving backwards.
        //    List<Point3D> backwardPathSegment = TraceStreamlineSegment(
        //        startPointOnMesh,
        //        stressVectors,
        //        effectiveStepLength,
        //        maximumTraceSteps,
        //        false); // false for negative direction

        //    // Combine the two segments into a single, continuous polyline.
        //    // The backward path is reversed (excluding the shared start point),
        //    // the common start point is added once, and then the forward path (excluding its start point).
        //    List<Point3D> completeStreamlinePoints = new List<Point3D>();

        //    // Add points from the backward trace, in reverse order, skipping the very first point
        //    // as it's the common start point and will be added explicitly.
        //    if (backwardPathSegment.Count > 1)
        //    {
        //        completeStreamlinePoints.AddRange(backwardPathSegment.Skip(1).Reverse());
        //    }

        //    // Add the shared starting point of the streamline.
        //    completeStreamlinePoints.Add(startPointOnMesh);

        //    // Add points from the forward trace, skipping the very first point.
        //    if (forwardPathSegment.Count > 1)
        //    {
        //        completeStreamlinePoints.AddRange(forwardPathSegment.Skip(1));
        //    }

        //    // Create and return the final polyline curve using the Synera geometry kernel.
        //    // A small tolerance (1E-05) is used for constructing the polyline.
        //    var polyline = Application.Current.KernelManager.Get<IGeometryKernel>().CreatePolyline(completeStreamlinePoints, 1E-05);
        //    return polyline;
        //}

        //private List<Point3D> TraceStreamlineSegment(
        //    Point3D segmentStartPoint,
        //    IDictionary<int, Vector3D> stressFieldVectors,
        //    double currentStepLength,
        //    int maxIterations,
        //    bool isMovingForward)
        //{
        //    // Initialize the list of points for this segment with the starting point.
        //    List<Point3D> segmentPoints = new List<Point3D> { segmentStartPoint };

        //    // The point from which the next step will be calculated.
        //    Point3D currentTracePoint = segmentStartPoint;

        //    // Determine the multiplier for the stress vector based on the desired tracing direction.
        //    int directionFactor = isMovingForward ? 1 : -1;

        //    // Loop for a maximum number of iterations to trace the streamline segment.
        //    for (int stepCounter = 0; stepCounter < maxIterations; ++stepCounter)
        //    {
        //        MeshParameter meshHitInfo;
        //        // Find the closest point on the mesh to the current tracing point, and get mesh component info.
        //        Point3D pointOnMeshFromCurrent = Mesh.ClosestPointOnMesh(currentTracePoint, out meshHitInfo);

        //        // Determine the original FEA element ID associated with the current mesh location.
        //        int currentElementId = GetElementIdFromMeshParameter(pointOnMeshFromCurrent, meshHitInfo);

        //        // If the element ID is invalid or no stress vector is found for it, stop tracing.
        //        if (currentElementId == -1 || !stressFieldVectors.ContainsKey(currentElementId))
        //        {
        //            break;
        //        }

        //        // Retrieve the base stress vector for the current element.
        //        Vector3D baseStressVector = stressFieldVectors[currentElementId];
        //        // Apply the direction factor to get the actual stepping vector.
        //        Vector3D effectiveStepVector = directionFactor * baseStressVector;

        //        // Ensure the direction of the current step is consistent with the previous segment.
        //        // This prevents jagged lines if the stress vector's orientation flips across element boundaries.
        //        if (segmentPoints.Count >= 2)
        //        {
        //            // Calculate the vector of the previously added segment.
        //            //Vector3D previousSegmentDirection = segmentPoints[segmentPoints.Count - 1].Subtract(segmentPoints[segmentPoints.Count - 2]);
        //            Vector3D previousSegmentDirection = Point3D.Subtract(segmentPoints[segmentPoints.Count - 1], segmentPoints[segmentPoints.Count - 2]);
        //            // If the dot product is negative, the vectors are pointing in opposite general directions, so reverse.
        //            if (Vector3D.DotProduct(effectiveStepVector, previousSegmentDirection) < 0.0)
        //            {
        //                effectiveStepVector.Reverse(); // Reverse the current step vector.
        //            }
        //        }

        //        // Calculate the next tentative point by moving along the effective stress vector.
        //        Point3D nextTentativeLocation = currentTracePoint + currentStepLength * effectiveStepVector;

        //        // Project this tentative point back onto the mesh to ensure the streamline stays on the surface.
        //        MeshParameter nextMeshHitInfo;
        //        Point3D nextPointOnMesh = Mesh.ClosestPointOnMesh(nextTentativeLocation, out nextMeshHitInfo);

        //        // Check if the streamline has hit a "naked" (boundary) edge of the mesh.
        //        // If so, add this point and terminate the segment tracing.
        //        if (nextMeshHitInfo.MeshComponent.Type == MeshComponentType.Edge && OpenMeshEdges.Contains(nextMeshHitInfo.MeshComponent.Index))
        //        {
        //            segmentPoints.Add(nextPointOnMesh);
        //            break;
        //        }

        //        // Check for stagnation: if the new point is effectively the same as the current point,
        //        // it means the streamline is stuck (e.g., at an extremum, or due to precision). Terminate.
        //        if (nextPointOnMesh.AlmostEquals(currentTracePoint, 1E-05)) // Using a small tolerance for floating-point comparison
        //        {
        //            break;
        //        }

        //        // Add the newly found point on the mesh to the segment.
        //        segmentPoints.Add(nextPointOnMesh);
        //        // Update the current tracing point for the next iteration.
        //        currentTracePoint = nextPointOnMesh;
        //    }

        //    return segmentPoints;
        //}

        //private int GetElementIdFromMeshParameter(Point3D pointOnMeshLocation, MeshParameter meshComponentInfo)
        //{
        //    int[] potentialFaceIndices; // Indices of mesh faces that are associated with the hit mesh component.

        //    // Determine the relevant face indices based on the type of mesh component hit.
        //    switch (meshComponentInfo.MeshComponent.Type)
        //    {
        //        case MeshComponentType.Vertex:
        //            potentialFaceIndices = Mesh.GetFaceIndicesForVertex(meshComponentInfo.MeshComponent.Index);
        //            break;
        //        case MeshComponentType.Edge:
        //            potentialFaceIndices = Mesh.GetFaceIndicesForEdge(meshComponentInfo.MeshComponent.Index);
        //            break;
        //        case MeshComponentType.Face:
        //            // If the hit is directly on a face, that's the only relevant face.
        //            potentialFaceIndices = new int[] { meshComponentInfo.MeshComponent.Index };
        //            break;
        //        case MeshComponentType.None:
        //        default:
        //            // If no valid mesh component is found, we cannot determine an element.
        //            return -1;
        //    }


        //    // If no associated faces were found for the mesh component, return -1.
        //    if (potentialFaceIndices == null || potentialFaceIndices.Length == 0)
        //    {
        //        return -1;
        //    }

        //    int determinedMeshFaceIndex;
        //    if (potentialFaceIndices.Length == 1)
        //    {
        //        // If only one face is associated, that's the one we'll use.
        //        determinedMeshFaceIndex = potentialFaceIndices[0];
        //    }
        //    else
        //    {
        //        // If multiple faces are associated (e.g., at a shared vertex or edge),
        //        // select the face whose center is geometrically closest to the hit point on the mesh.
        //        determinedMeshFaceIndex = potentialFaceIndices
        //            .OrderBy(faceIdx => pointOnMeshLocation.SquaredDistanceTo(Mesh.GetFaceCenter(faceIdx)))
        //            .First();
        //    }

        //    // Map the selected mesh face index back to its original FEA element ID using the stored mapping.
        //    if (determinedMeshFaceIndex >= 0 && determinedMeshFaceIndex < FaceToElementIndices.Count)
        //    {
        //        return FaceToElementIndices[determinedMeshFaceIndex];
        //    }
        //    return -1; // Indicate that no valid FEA element could be mapped.
        //}
    }
}
