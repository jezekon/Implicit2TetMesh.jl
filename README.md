# Implicit2TetMesh.jl

Implicit2TetMesh is an experimental Julia package for generating high-quality tetrahedral meshes from implicit geometries defined by Signed Distance Functions (SDFs). Implementation details are provided in the documentation below. For practical usage examples, see test/Examples.

<div style="display: flex; justify-content: center; align-items: center; gap: 10px;">
  <img src="doc/beam.png" style="height: 270px; max-width: 50%;" alt="Original beam geometry" />
  <img src="doc/beam_cut.png" style="height: 270px; max-width: 50%;" alt="Sliced beam tetrahedral mesh" />
</div>

## Features
- **Robust Meshing**: High-quality tetrahedral mesh generation from implicit geometries using A15 (body-centered cubic) and Schlafli (orthoscheme) discretizations
- **Isosurface Refinement**: Advanced boundary processing with experimental NZZZ case handling for thin features relative to characteristic element size
- **Volume Correction**: Precise volume preservation using automatic surface node adjustment
- **Geometric Constraints**: Bounded plane definitions for selective node alignment
- **Mesh Operations**: Slicing, isolated component removal, inverted element fixing, and VTU export with mesh quality metrics

## Function Input
Main function for generating tetrahedral meshes from SDF data:

```julia
generate_tetrahedral_mesh(grid_file, sdf_file, output_prefix; options=MeshGenerationOptions())
```

### Parameters:
- `grid_file::String`: Path to the JLD2 file containing the grid data
- `sdf_file::String`: Path to the JLD2 file containing the SDF values
- `output_prefix::String`: Prefix for output files (default: "output")
- `options::MeshGenerationOptions`: Configuration options (optional)

### Return Value:
- `mesh::BlockMesh`: The generated tetrahedral mesh
- `vtu`: Mesh visualization files in Paraview format

## MeshGenerationOptions
The `MeshGenerationOptions` struct allows for customization of the mesh generation process:

```julia
MeshGenerationOptions(;
    scheme="A15",                          # Discretization scheme: "A15" or "Schlafli"
    warp_param=0.3,                        # Parameter controlling warping behavior
    plane_definitions=nothing,             # Optional cutting planes
    quality_export=false,                  # Whether to export mesh with quality metrics
    split_elements=true                    # Whether to split elements along the isosurface
)
```

### Options:
#### scheme::String
- Discretization scheme to use. Valid values:
  - `"A15"` (default): Higher quality tetrahedral elements
  - `"Schlafli"`: Alternative scheme with different element patterns

#### warp_param::Float64
- Parameter controlling the warping behavior for cutting planes
- Must be non-negative (≥ 0)
- Default: `0.3`

#### plane_definitions::Union{Vector{PlaneDefinition}, Nothing}
- Optional cutting planes for mesh modification
- Use `PlaneDefinition` objects to define planes with different shapes
- Default: `nothing` (no cutting planes)

#### quality_export::Bool
- Controls whether detailed quality metrics are included in output files
- Default: `false`

#### split_elements::Bool
- Controls element handling at the isosurface boundary:
  - `true` (default): Split elements along the isosurface for higher accuracy
  - `false`: Move nodes to the isosurface without splitting elements

## Plane Definitions
For creating cutting planes to modify the mesh:

```julia
PlaneDefinition(normal, point, shape)
```

### Parameters:
- `normal::Vector{Float64}`: Normal vector of the plane
- `point::Vector{Float64}`: Point on the plane
- `shape::PlaneShape`: Shape of the bounded plane (Rectangle, Square, Circle, or Ellipse)

### Plane Shapes:
```julia
Rectangle(width, height)    # Rectangular bounded plane
Square(size)                # Square bounded plane (special case of Rectangle)
Circle(radius)              # Circular bounded plane
Ellipse(a, b)               # Elliptical bounded plane
```

## Example Usage

```julia
using Implicit2TetMesh

# Basic usage with default options
mesh = generate_tetrahedral_mesh(
    "path/to/grid_data.jld2",
    "path/to/sdf_data.jld2",
    "beam"
)
```
## Advanced Example Usage

```julia
# Advanced usage with custom options and cutting planes
plane_definitions = [
    PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
    PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0))
]

options = MeshGenerationOptions(
    scheme = "A15",
    warp_param = 0.5,
    plane_definitions = plane_definitions,
    quality_export = true,
    split_elements = true
)

mesh = generate_tetrahedral_mesh(
    "path/to/grid_data.jld2",
    "path/to/sdf_data.jld2",
    "beam_cut",
    options=options
)
slice_mesh_with_plane!(mesh, "x", 0.5, -1) # sliced by yz plane in half, -1 (normal) -> delete every element before plane
export_mesh_vtu(mesh, "sliced.vtu")
```

## Quality Assessment

The package provides comprehensive mesh quality assessment through:

```julia
assess_mesh_quality(mesh, output_prefix)
```

This function evaluates:
- Geometric quality metrics (aspect ratio, dihedral angles, shape quality)
- Element validity (inverted elements)
- Boundary integrity
- Element orientation consistency
- Element intersections
- Dihedral angle distribution with histogram visualization

For quick volume checks:

```julia
TetMesh_volumes(mesh)
```

## TODO List
- [ ] Parallel processing for large meshes
