# Implicit2TetMesh.jl

A Julia package for generating high-quality tetrahedral meshes from implicit geometries defined by Signed Distance Functions (SDFs). It transforms smooth implicit surfaces into finite element tetrahedral mesh while preserving geometric features and ensuring element quality. The package includes robust mesh optimization, quality assessment, and geometric modification capabilities.

<div style="display: flex; justify-content: center; align-items: center; gap: 10px;">
  <img src="doc/beam.png" style="height: 250px; max-width: 50%;" alt="Original beam geometry" />
  <img src="doc/beam_cut.png" style="height: 250px; max-width: 50%;" alt="Sliced beam tetrahedral mesh" />
</div>

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
    optimize=true,                         # Whether to perform mesh optimization
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

#### optimize::Bool
- Determines whether mesh optimization is performed
- Default: `true`

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
    optimize = true,
    split_elements = true
)

mesh = generate_tetrahedral_mesh(
    "path/to/grid_data.jld2",
    "path/to/sdf_data.jld2",
    "beam_cut",
    options=options
)

# Without mesh optimization
mesh = generate_tetrahedral_mesh(
    "path/to/grid_data.jld2",
    "path/to/sdf_data.jld2",
    "beam_no_opt",
    options=MeshGenerationOptions(optimize = false)
)
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

## Features

- **Robust Meshing**: Creates high-quality tetrahedral meshes from implicit surfaces using Implicit Domain Stuffing algorithm
- **Mesh Optimization**: Includes Laplacian smoothing and adaptive quality-based optimization
- **Quality Control**: Quality metrics and visualization tools
- **Geometric Modifications**: Support for cutting planes and boundary refinement
- **Multiple Discretization Schemes**: A15 and Schlafli schemes for different element patterns
- **Export Capabilities**: VTK-based visualization with detailed quality metrics

## TODO List
- [ ] Parallel processing for large meshes
- [ ] Further optimization of element quality metrics
