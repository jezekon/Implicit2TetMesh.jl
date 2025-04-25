"""
    MeshGenerationOptions

Configuration options for tetrahedral mesh generation.

# Fields
- `scheme::String`: Discretization scheme, either "A15" or "Schlafli" (default: "A15")
- `warp_param::Float64`: Parameter controlling warping behavior (default: 0.3)
- `plane_definitions::Union{Vector{PlaneDefinition}, Nothing}`: Optional cutting planes (default: nothing)
- `optimize::Bool`: Whether to perform mesh optimization (default: true)
- `split_elements::Bool`: Whether to split elements along the isosurface (true) or just move nodes to the isosurface (false) (default: true)
"""
struct MeshGenerationOptions
    scheme::String
    warp_param::Float64
    plane_definitions::Union{Vector{PlaneDefinition}, Nothing}
    optimize::Bool
    split_elements::Bool
    
    # Inner constructor with defaults
    function MeshGenerationOptions(;
        scheme::String = "A15",
        warp_param::Float64 = 0.3,
        plane_definitions::Union{Vector{PlaneDefinition}, Nothing} = nothing,
        optimize::Bool = true,
        split_elements::Bool = true
    )
        # Validate scheme selection
        if !(scheme in ["A15", "Schlafli"])
            error("Invalid scheme: $scheme. Must be 'A15' or 'Schlafli'")
        end

        # Validate warp_param to be non-negative (kladný nebo nulový)
        if warp_param < 0.0
            error("Invalid warp_param: $warp_param. Must be non-negative (>= 0).")
        end
        
        new(scheme, warp_param, plane_definitions, optimize, split_elements)
    end
end

"""
    generate_tetrahedral_mesh(grid_file::String, sdf_file::String, output_prefix::String = "output"; 
                             options::MeshGenerationOptions = MeshGenerationOptions())

Generate a tetrahedral mesh from implicit surface definition (SDF) data.

# Arguments
- `grid_file::String`: Path to the JLD2 file containing the grid data
- `sdf_file::String`: Path to the JLD2 file containing the SDF values
- `output_prefix::String`: Prefix for output files (default: "output")
- `options::MeshGenerationOptions`: Configuration options for mesh generation

# Returns
- `mesh::BlockMesh`: The generated tetrahedral mesh

# Example
```julia
# Basic usage with default options
mesh = generate_tetrahedral_mesh("grid_data.jld2", "sdf_data.jld2", "beam")

# Advanced usage with custom options
options = MeshGenerationOptions(
    scheme = "Schlafli",
    warp_param = 0.5,
    plane_definitions = [
        PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
        PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0))
    ]
)
mesh = generate_tetrahedral_mesh("grid_data.jld2", "sdf_data.jld2", "beam_cut", options=options)
"""

function generate_tetrahedral_mesh(grid_file::String, sdf_file::String, output_prefix::String = "output";
    options::MeshGenerationOptions = MeshGenerationOptions())
    # Step 1: Load data from JLD2 files
    @info "Loading data from $grid_file and $sdf_file..."
    @load grid_file fine_grid
    @load sdf_file fine_sdf
    
    # Step 2: Create the BlockMesh from the loaded data
    @info "Creating BlockMesh structure..."
    mesh = BlockMesh(fine_sdf, fine_grid)
    
    # Step 3: Generate mesh with the chosen discretization scheme
    @info "Generating mesh with $(options.scheme) scheme..."
    generate_mesh!(mesh, options.scheme)
    
    # Step 4: Warp nodes to the isocontour (zero level set of SDF)
    @info "Warping nodes to isocontour..."
    warp!(mesh)
    
    # Step 5: Update mesh connectivity (cleanup nodes, merge duplicates)
    @info "Updating mesh connectivity..."
    update_connectivity!(mesh)
    
    # Step 6: Process the isosurface boundary using the selected method
    if options.split_elements
        @info "Slicing ambiguous tetrahedra (with element splitting)..."
        slice_ambiguous_tetrahedra!(mesh)
    else
        @info "Adjusting nodes to isosurface (without element splitting)..."
        adjust_nodes_to_isosurface!(mesh)
    end
    
    # Step 7: Update mesh connectivity again after isosurface processing
    @info "Updating mesh connectivity..."
    update_connectivity!(mesh)
    
    # Step 8: Compute and display tetrahedra volume information
    @info "Computing mesh volumes..."
    TetMesh_volumes(mesh)
    
    # Step 9: Optimize the mesh if requested
    if options.optimize
        @info "Optimizing mesh..."
        optimize_mesh!(mesh)
    end
    
    # Step 10: Export the initial mesh to VTK format
    output_file = "$(output_prefix)_TriMesh-$(options.scheme).vtu"
    @info "Exporting mesh to $output_file..."
    export_mesh_vtk(mesh, output_file)
    
    # Step 11: Apply cutting planes if defined
    if options.plane_definitions !== nothing && options.warp_param !== 0.
        @info "Applying cutting planes with warp_param = $(options.warp_param)..."
        warp_mesh_by_planes_sdf!(mesh, options.plane_definitions, options.warp_param)
        
        # Update mesh connectivity after applying cutting planes
        @info "Updating mesh connectivity..."
        update_connectivity!(mesh)
        
        # Export the cut mesh to VTK
        cut_output_file = "$(output_prefix)_TriMesh-$(options.scheme)_cut.vtu"
        @info "Exporting cut mesh to $cut_output_file..."
        export_mesh_vtk(mesh, cut_output_file)
    end
    
    return mesh
end

"""
## Usage Examples
# Example 1: Basic usage with default options
mesh = generate_tetrahedral_mesh(
    "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2",
    "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2",
    "cantilever_beam_interp"
)

# Example 2: Using Schlafli scheme with custom cutting planes
plane_definitions = [
    PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
    PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0))
]

mesh = generate_tetrahedral_mesh(
    "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2",
    "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2",
    "cantilever_beam_interp_cut",
    options=MeshGenerationOptions(
        scheme = "Schlafli",
        warp_param = 0.5,
        plane_definitions = plane_definitions
    )
)

# Example 3: Without mesh optimization
mesh = generate_tetrahedral_mesh(
    "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2",
    "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2",
    "cantilever_beam_interp_nopt",
    options=MeshGenerationOptions(optimize = false)
)
"""
