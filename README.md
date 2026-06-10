# Implicit2TetMesh.jl

Implicit2TetMesh is an experimental Julia package for generating high-quality tetrahedral meshes from implicit geometries defined by Signed Distance Functions (SDFs), inspired by [isosurface stuffing algorithm](https://dl.acm.org/doi/10.1145/1276377.1276448). Implementation details are provided in the documentation below. For practical usage examples, see [`test/Examples/`](test/Examples/).

<!-- <div style="display: flex; justify-content: center; align-items: center; gap: 10px;">
  <img src="doc/beam.png" style="height: 270px; max-width: 50%;" alt="Original beam geometry" />
  <img src="doc/beam_cut.png" style="height: 270px; max-width: 50%;" alt="Sliced beam tetrahedral mesh" />
</div> -->
<p align="center">
  <img src="doc/gripper-half_sdf_half_tet.png" width="80%" alt="Raw topology optimization result" />
</p>

## Features
- **Robust Meshing**: High-quality tetrahedral mesh generation from implicit geometries on an A15 (body-centered cubic) acute lattice
- **Edge-based Warping**: Lattice vertices near the surface are snapped onto the linear cut points of sign-crossing edges (isosurface-stuffing warp), which collapses the sliver tetrahedra
- **Boundary Trimming**: Tetrahedra crossing the surface are trimmed back to the interior with quartet-style `trim_spikes` stencils; a consistent quad-split diagonal keeps the trimmed boundary crack-free. Surface-only ("quadruple-zero") tetrahedra are resolved with Labelle's §3.4 heuristic, which removes bubbles at the trimming stage instead of relying on a global cleanup pass
- **Geometric Constraints**: Bounded plane definitions for selective node alignment
- **Mesh Operations**: Slicing, isolated component removal, inverted element fixing, and VTU export with mesh quality metrics

## Installation

**Requirements:** Julia LTS (1.10.10)
```julia
# From Julia REPL, press ] to enter package mode
pkg> add https://github.com/jezekon/Implicit2TetMesh.jl
```
or
```
git clone https://github.com/jezekon/Implicit2TetMesh.jl
```
## Main Function
Generate tetrahedral meshes from SDF data:
```julia
generate_tetrahedral_mesh(grid_file, sdf_file, output_prefix; options=MeshGenerationOptions())
```
#### Parameters:
- `grid_file::String`: Path to the JLD2 file containing the grid data
- `sdf_file::String`: Path to the JLD2 file containing the SDF values
- `output_prefix::String`: Prefix for output files (default: "output")
- `options::MeshGenerationOptions`: Configuration options (optional)
#### SDF Convention:
The package uses the standard convention `phi < 0 = inside`, `phi > 0 = outside`, `phi = 0 = on the surface`, matching the isosurface-stuffing reference implementations. Input files that store the opposite sign (positive = inside) are negated automatically when the `BlockMesh` is constructed; the data files on disk are never modified.
#### Surface Warping:
Before the boundary is sliced, lattice vertices close to the surface are snapped onto the **linear cut points** of their sign-crossing edges, following the isosurface-stuffing warp (Labelle 2007, §3.2). For every tetrahedron edge whose endpoints have opposite SDF signs, the surface crosses the edge at `X_i + α·(X_j − X_i)` with `α = φ_i / (φ_i − φ_j)`; if that crossing lies within a fraction `threshold` (default `0.3`) of an endpoint, the endpoint is moved onto it and marked as lying on the surface. Each vertex is warped to its closest qualifying cut point, and all displacements are computed from the original geometry and applied at once, so the result is independent of vertex order. Snapping to the *linear* cut point (rather than projecting each node to the true isosurface) is what keeps the surrounding tetrahedra well shaped and removes the slivers.
#### Surface Tetrahedron Handling:
After trimming, some tetrahedra end up with **all four vertices on the surface** (every node was warped onto a cut point). These "quadruple-zero" tetrahedra are ambiguous — it is not clear whether each one lies inside or outside the geometry — so they are resolved with the principled heuristic of Labelle 2007 (§3.4) rather than a single centroid test:

1. A candidate that is **inverted**, or whose **dihedral angles** fall outside `[10°, 140°]`, is discarded — it is too flat to improve surface fidelity.
2. Of the well-shaped survivors, the number of faces that **adjoin the interior mesh** decides the outcome: all four faces adjoin → **retain** (the tetrahedron fills a pocket in the boundary); no face adjoins → **discard** (an isolated "bubble"); otherwise the **SDF sign at the centroid** decides (inside → keep).

Because bubbles are removed here, the global `remove_isolated_components!` pass that follows is now only a safety net (it keeps the largest connected component for geometries that split into disconnected pieces). Where four warped vertices lie on a sub-lattice-thickness thin feature, the boundary can still self-touch (two surface sheets meeting along an edge); such pinches reflect the uniform lattice resolution and are the target of the planned adaptive refinement, not cracks — the trimmed boundary stays free of open (hole) edges.
#### Volume Accuracy:
The mesh volume is an honest discretization of the zero isosurface. Because the boundary is approximated by flat triangles, the meshed volume is slightly smaller than the reference SDF volume wherever the isosurface is curved (about 2.7 % on the beam, 0.4 % on the gripper); this gap shrinks as the lattice is refined, not by moving nodes after meshing.

The package therefore applies **no** post-hoc volume correction. Displacing surface nodes onto a single global SDF level after meshing is equivalent to meshing the `phi = c` isocontour, but done crudely — it degrades the surface fidelity the warp establishes and can produce spiked or inverted elements (the quartet / Labelle reference algorithms have no such step). If an exact target volume is ever required (for example a volume fraction carried over from topology optimization), the recommended approach is to choose the iso-level offset `c` by bisection **before** meshing and run the normal pipeline on the shifted field `phi − c`. That keeps the mesh robust (no spikes or inversions) and preserves the full surface-fidelity guarantee.
#### Return Value:
- `mesh::BlockMesh`: The generated tetrahedral mesh
- **Output files**: `.vtu` mesh visualization files for Paraview

### MeshGenerationOptions

Configure the mesh generation process with the following options:

```julia
MeshGenerationOptions(;
    scheme::String = "A15",                           # Discretization scheme (only "A15" is supported)
    warp_param::Float64 = 0.3,                        # Warping intensity for plane alignment (0.0 = disabled)
    plane_definitions::Union{Vector{PlaneDefinition}, Nothing} = nothing,  # Cutting planes for BC application
    quality_export::Bool = false                      # Export detailed quality metrics
)
```
#### Option Details

- **scheme**: `"A15"` — body-centered cubic acute lattice (the only supported scheme)
- **warp_param**: Controls how strongly nodes are attracted to cutting planes (0.0-1.0 range recommended)
- **plane_definitions**: Vector of `PlaneDefinition` objects for boundary plane constraints
- **quality_export**: When `true`, exports additional quality metrics (Jacobian determinants, dihedral angles, volume ratios)

### Example Usage
```julia
using Implicit2TetMesh

# Basic usage with default options
mesh = generate_tetrahedral_mesh(
    "path/to/grid_data.jld2",
    "path/to/sdf_data.jld2",
    "beam"
)
```
### Advanced Usage Examples
For complete examples with detailed documentation, see [`test/Examples/`](test/Examples/):
```julia
# Run beam example
julia --project=. test/Examples/beam.jl

julia --project=. test/Examples/gripper.jl
```
___
## TODO List
- [x] Principled surface-tetrahedron (quadruple-zero) handling per Labelle §3.4
- [ ] Optional adaptive refinement for thin features (sub-lattice-thickness walls)
- [ ] Performance optimizations for large meshes

## Acknowledgments
This package is an experimental tool. Please validate results carefully.
