# AidaTool — Stress-Guided Mesh Deformation for Synera

AidaTool is a C# addin for [Synera](https://www.synera.io/) that deforms finite element shell meshes so their topology conforms to the internal stress state of a structure. The goal is a mesh whose elements are elongated and aligned along the principal stress trajectories — the directions in which the structure is most strongly loaded.

This concept is inspired by **graphic statics**, where the geometry of a structure (the form diagram) is reciprocal to the distribution of internal forces (the force diagram).

---

## Nodes

| Node | File | Mesh Type | Update Scheme | Parameters |
|------|------|-----------|---------------|------------|
| `FieldCreation` | `FieldCreation.cs` | Any | Gauss-Seidel | Iterations, Alpha |
| `QuadMeshFD` | `QuadMeshFD.cs` | Quad | Gauss-Seidel | Iterations, Alpha, Beta |
| `TriMeshFD` | `TriMeshFD.cs` | Tri | Jacobi + Lambda | Iterations, Alpha, Beta, Lambda |
| `MeshFD` | `MeshFD.cs` | Quad / Tri / Mixed | Auto-detected | Iterations, Alpha, Beta, Lambda |
| `MeshFD3D` | `MeshFD3D.cs` | Quad / Tri / Mixed | Auto-detected | Iterations, Alpha, Beta, Lambda + Target Surface |

---

## Background: Laplacian Smoothing

Standard **Laplacian smoothing** repositions each interior node by moving it to the geometric centroid of its neighbours:

```
x_i  ←  (1 / N) * Σ x_j      for all neighbours j of node i
```

Applied iteratively, this distributes nodes evenly, reducing element distortion. It has no awareness of the stress field — every neighbour is weighted equally regardless of load level or direction.

---

## Core Algorithm (all nodes except FieldCreation)

Each relaxation iteration applies three steps per interior node.

### Step 1 — Von Mises Weighted Centroid

Neighbours in high-stress regions attract the current node more strongly:

```
weight_j  =  1 + α * σ̄_j

centroid  =  Σ (weight_j * x_j) / Σ weight_j

δ  =  centroid − x_i
```

Where `σ̄_j` is the von Mises stress at neighbour `j`, normalized globally to `[0, 1]`.

**Effect:** Interior nodes migrate toward high-stress regions. The mesh becomes stress-aware in terms of *where* stresses are high.

### Step 2 — Decompose Displacement into Principal Stress Components

The displacement vector `δ` is decomposed into three orthogonal components:

```
δ_major  =  (δ · ê_major) * ê_major       [along major principal stress]
δ_minor  =  (δ · ê_minor) * ê_minor       [along minor principal stress]
δ_rest   =  δ − δ_major − δ_minor         [residual / out-of-plane]
```

Where `ê_major` and `ê_minor` are the local principal stress directions, averaged over all elements connected to the node.

> **Why decompose the displacement vector instead of weighting neighbours directly?**
> If neighbours lie on both sides of the stress axis, their directional contributions cancel exactly — the node does not move preferentially. Decomposing the *net* displacement vector avoids this cancellation entirely.

### Step 3 — Amplify Stress-Aligned Components and Preserve Step Length

The stress-aligned components are amplified by `β` and the local dominance of each direction:

```
δ_anisotropic  =  δ_rest
               +  (1 + β * n_major) * δ_major
               +  (1 + β * n_minor) * δ_minor
```

Where `n_major` and `n_minor` are **per-node normalized** principal stress magnitudes:

```
localMax  =  max(|σ_major|, |σ_minor|)  at node i
n_major   =  |σ_major| / localMax
n_minor   =  |σ_minor| / localMax
```

Per-node normalization ensures the dominant principal direction always contributes with `n = 1.0`, regardless of how small the absolute stress is compared to other parts of the model.

**Step length is always preserved.** The anisotropic vector is rescaled to the original displacement length before being applied:

```
x_i  ←  x_i  +  |δ| * (δ_anisotropic / |δ_anisotropic|)
```

`β` only steers the *direction* of movement — it never changes the step size. Mesh inversion is not possible regardless of how large `β` is set.

---

## Update Schemes: Gauss-Seidel vs. Jacobi

### Gauss-Seidel (QuadMeshFD, MeshFD quad/mixed path)

Each node is updated **in-place** as the loop progresses. Subsequent nodes immediately see the already-updated positions of earlier nodes. This converges quickly and is stable for quad meshes, which typically have ~4 neighbours at 90° angles.

### Jacobi (TriMeshFD, MeshFD tri path)

Before each iteration, a **frozen snapshot** of all current positions is taken. Each node is computed using only snapshot values and all updates are written simultaneously at the end of the pass. This prevents chaotic mid-iteration propagation through the dense, tightly angled connectivity of triangular meshes (~6 neighbours at 60°).

### Lambda — Damping for Triangular Meshes

Even with Jacobi updates, triangular meshes can develop sliver elements if per-iteration displacements are too large. Lambda `(0, 1]` scales down each displacement step:

```
x_i  ←  x_i  +  λ * δ_anisotropic_rescaled
```

`λ = 1.0` (default) applies the full step. Lower values slow convergence but prevent element collapse.

Lambda is only active on the triangular path. For quad and mixed meshes it is forced to `1.0`.

---

## 3D Surface Projection (MeshFD3D only)

When smoothing a mesh that lies on a **curved 3D surface**, the weighted average of neighbouring positions lies *below* the surface (inside the curvature) — not on it. This drift compounds over iterations, pulling nodes off the geometry.

MeshFD3D corrects this by adding a projection step **after every displacement**:

```
candidatePosition  =  x_i + λ * δ_anisotropic_rescaled

x_i  ←  targetSurface.ClosestPoint(candidatePosition)
```

The candidate position is snapped to the nearest point on the reference surface before being committed. This happens every step, so drift never accumulates.

---

## Parameters

| Parameter | Type | Range | Default | Effect |
|-----------|------|-------|---------|--------|
| **Iterations** | Integer | ≥ 1 | — | Number of relaxation passes. Typical range: 5–50. |
| **Alpha** | Integer | ≥ 0 | — | Von Mises stress influence. `α = 0` gives uniform Laplacian smoothing. Typical range: 1–10. |
| **Beta** | Double | ≥ 0 | — | Principal stress direction bias. `β = 0` disables directional anisotropy. Typical range: 0.5–5. |
| **Lambda** | Double | (0, 1] | 1.0 | Damping for triangular meshes. Inactive (forced to 1.0) for quad/mixed meshes. |

---

## Inputs and Outputs

### FieldCreation

| Port | Direction | Type | Description |
|------|-----------|------|-------------|
| Input model | In | `IModel` | Solved FEA model with shell elements and stress results. |
| Iterations | In | `SyneraInt` | Number of relaxation iterations. |
| Alpha | In | `SyneraInt` | Von Mises stress influence strength. |
| Von-Mises Stress | Out | `List<double>` | Von Mises stress per node. |
| Vertices | Out | `List<Point3D>` | Original node positions. |
| Deformed Vertices | Out | `List<Point3D>` | Node positions after deformation. |
| Mesh | Out | `IMesh` | Deformed mesh. |

### QuadMeshFD / TriMeshFD

Same as FieldCreation plus:

| Port | Direction | Type | Description |
|------|-----------|------|-------------|
| Beta | In | `SyneraDouble` | Principal stress direction bias. |
| Lambda | In | `SyneraDouble` | Damping factor — TriMeshFD only. |
| Max Principle Stress | Out | `List<Vector3D>` | Major principal stress vectors per element. |
| Min Principle Stress | Out | `List<Vector3D>` | Minor principal stress vectors per element. |

### MeshFD

Same as QuadMeshFD plus Lambda input. Auto-detects mesh type and selects update scheme accordingly.

### MeshFD3D

Same as MeshFD plus:

| Port | Direction | Type | Description |
|------|-----------|------|-------------|
| Target surface | In | `IMesh` | 3D reference surface onto which all nodes are projected after each displacement step. |

---

## Mesh Type Detection (MeshFD and MeshFD3D)

| Mesh | Detection | Update Scheme | Lambda |
|------|-----------|---------------|--------|
| Pure triangular | All elements have 3 corner nodes | Jacobi | User-defined |
| Pure quad | All elements have 4 corner nodes | Gauss-Seidel | Forced to 1.0 |
| Mixed | Both 3- and 4-corner elements present | Gauss-Seidel (quad-dominant) | Forced to 1.0 |

---

## Implementation Notes

### Stress Data Pipeline

1. Element-level principal stress vectors are computed by `StressStreamlineBuilder.GetPrincipalStressValues()`. Vectors are scaled by eigenvalue — length encodes intensity, direction encodes axis.
2. Both von Mises scalar and principal vectors are averaged over all elements connected to each node (corner nodes only; mid-side nodes excluded).
3. Von Mises stress is computed from the full 6-component tensor `[Sxx, Syy, Szz, Sxy, Syz, Sxz]`:

```
σ_VM  =  √[ 0.5 * ((Sxx−Syy)² + (Syy−Szz)² + (Szz−Sxx)²  +  6(Sxy² + Syz² + Sxz²)) ]
```

### Boundary Conditions

Nodes on naked (free) edges are identified at the start and kept fixed throughout all iterations. Only interior nodes participate in the relaxation. This preserves the outer shape and footprint of the mesh.

### Normalization

- **Von Mises stress** — normalized globally (min-max across all nodes). `α` is scale-independent regardless of unit system or load magnitude.
- **Principal stress magnitudes** — normalized per-node using the larger of the two local magnitudes. Prevents a single stress peak from collapsing the directional contribution of all other nodes to near-zero.

---

## Node Comparison

| Feature | FieldCreation | QuadMeshFD | TriMeshFD | MeshFD | MeshFD3D |
|---------|:---:|:---:|:---:|:---:|:---:|
| Von Mises stress weighting (Alpha) | Yes | Yes | Yes | Yes | Yes |
| Principal stress direction bias (Beta) | No | Yes | Yes | Yes | Yes |
| Gauss-Seidel update | Yes | Yes | No | Auto | Auto |
| Jacobi update | No | No | Yes | Auto | Auto |
| Lambda damping | No | No | Yes | Yes | Yes |
| Mixed mesh support | No | No | No | Yes | Yes |
| 3D surface projection | No | No | No | No | Yes |

Setting `β = 0` in any of the FD nodes produces output identical to `FieldCreation`. Increasing `β` progressively elongates elements along the principal stress trajectories while keeping overall mesh density controlled by `α`.
