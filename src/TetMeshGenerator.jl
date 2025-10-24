"""
    MeshGenerationOptions

Configuration options for tetrahedral mesh generation.

# Fields
- `scheme::String`: Discretization scheme ("A15" or "Schlafli", default: "A15")
- `warp_param::Float64`: Warping intensity for surface nodes (default: 0.3)
- `plane_definitions::Union{Vector{PlaneDefinition}, Nothing}`: Cutting plane constraints (optional)
- `quality_export::Bool`: Export detailed quality metrics (default: false)
- `correct_volume::Bool`: Apply volume correction to match SDF reference (default: false)
- `experimental_nzzz::Bool`: Enable experimental NZZZ case processing (default: false)
"""
struct MeshGenerationOptions
    scheme::String
    warp_param::Float64
    plane_definitions::Union{Vector{PlaneDefinition},Nothing}
    quality_export::Bool
    correct_volume::Bool
    experimental_nzzz::Bool

    function MeshGenerationOptions(;
        scheme::String = "A15",
        warp_param::Float64 = 0.3,
        plane_definitions::Union{Vector{PlaneDefinition},Nothing} = nothing,
        quality_export::Bool = false,
        correct_volume::Bool = false,
        experimental_nzzz::Bool = false,
    )
        # Validate inputs
        scheme in ["A15", "Schlafli"] ||
            error("Invalid scheme: $scheme. Must be 'A15' or 'Schlafli'")
        warp_param >= 0.0 || error("Invalid warp_param: $warp_param. Must be non-negative.")

        new(
            scheme,
            warp_param,
            plane_definitions,
            quality_export,
            correct_volume,
            experimental_nzzz,
        )
    end
end

"""
    generate_tetrahedral_mesh(grid_file, sdf_file, output_prefix="output"; 
                              options=MeshGenerationOptions())

Generate a tetrahedral mesh from SDF (Signed Distance Function) data with optional 
plane constraints.

# Arguments
- `grid_file::String`: Path to JLD2 file with grid coordinates
- `sdf_file::String`: Path to JLD2 file with SDF values
- `output_prefix::String`: Output filename prefix (default: "output")
- `options::MeshGenerationOptions`: Generation configuration

# Returns
- `BlockMesh`: Generated tetrahedral mesh

# Example
```julia
# Basic usage
mesh = generate_tetrahedral_mesh("grid.jld2", "sdf.jld2", "beam")

# With cutting planes
planes = [PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0))]
options = MeshGenerationOptions(
    scheme = "Schlafli",
    warp_param = 0.5,
    plane_definitions = planes
)
mesh = generate_tetrahedral_mesh("grid.jld2", "sdf.jld2", "beam_cut", options=options)
```
"""
function generate_tetrahedral_mesh(
    grid_file::String,
    sdf_file::String,
    output_prefix::String = "output";
    options::MeshGenerationOptions = MeshGenerationOptions(),
)
    # Load SDF and grid data from JLD2 files
    @info "Loading data from $grid_file and $sdf_file..."
    @load grid_file fine_grid
    @load sdf_file fine_sdf

    # Initialize mesh structure from SDF data
    mesh = BlockMesh(fine_sdf, fine_grid)

    # Generate base tetrahedral mesh using selected discretization scheme
    generate_mesh!(mesh, options.scheme)

    # Warp nodes to isosurface (SDF = 0 level set)
    warp!(mesh, options.scheme)
    update_connectivity!(mesh)

    # Process isosurface boundary - remove exterior elements
    slice_ambiguous_tetrahedra!(mesh, options.scheme, options.experimental_nzzz)
    update_connectivity!(mesh)
    remove_inverted_elements!(mesh)

    # Remove disconnected/isolated mesh components
    remove_isolated_components!(mesh, keep_largest = true)
    update_connectivity!(mesh)

    # Display volume statistics
    @info "Computing mesh volumes..."
    TetMesh_volumes(mesh)

    # Optional: correct volume to match SDF reference
    if options.correct_volume
        correct_mesh_volume!(
            mesh,
            fine_sdf,
            fine_grid,
            options.scheme,
            plane_definitions = options.plane_definitions,
        )
    end

    # Export initial mesh to VTK format
    output_file = "$(output_prefix)_TriMesh-$(options.scheme).vtu"
    @info "Exporting mesh to $output_file..."
    export_mesh(mesh, output_file, options.quality_export)

    # Optional: apply cutting plane constraints
    if options.plane_definitions !== nothing && options.warp_param !== 0.0
        @info "Applying cutting planes with warp_param = $(options.warp_param)..."
        warp_mesh_by_planes_sdf!(mesh, options.plane_definitions, options.warp_param)

        # Export mesh with applied plane constraints
        cut_output_file = "$(output_prefix)_TriMesh-$(options.scheme)_cut.vtu"
        @info "Exporting cut mesh to $cut_output_file..."
        export_mesh(mesh, cut_output_file, options.quality_export)
    end

    return mesh
end

"""
    export_mesh(mesh, filename, quality_export)

Export tetrahedral mesh to VTK format with optional quality metrics.

# Arguments
- `mesh::BlockMesh`: Mesh to export
- `filename::String`: Output filename (with .vtu extension)
- `quality_export::Bool`: Include detailed quality metrics if true
"""
function export_mesh(mesh::BlockMesh, filename::String, quality_export::Bool)
    export_func = quality_export ? export_mesh_vtu_quality : export_mesh_vtu
    return export_func(mesh, filename)
end
