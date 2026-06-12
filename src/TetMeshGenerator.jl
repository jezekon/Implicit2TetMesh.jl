"""
    MeshGenerationOptions

Configuration options for tetrahedral mesh generation.

# Fields
- `scheme::String`: Discretization scheme (only "A15" is supported, default: "A15")
- `warp_param::Float64`: Warping intensity for surface nodes (default: 0.3)
- `plane_definitions::Union{Vector{PlaneDefinition}, Nothing}`: Cutting plane constraints (optional)
- `quality_export::Bool`: Export detailed quality metrics (default: false)
- `cut_points::Symbol`: How the surface cut point on a sign-crossing edge is located, in both
  the warp and the slicing stage (default: `:linear`).
    * `:linear` -- quartet's estimate from the two endpoint SDF values. Exact when the input
      is a true signed distance function (the field is then linear along lattice edges) and
      keeps bit-identity with the quartet reference implementation.
    * `:bisection` -- the true zero of the interpolated field along the edge (Labelle &
      Shewchuk 2007 §3.1). Use this when the input field is NOT distance-like (e.g. an
      RBF-smoothed SDF): with `:linear` such fields get surface vertices misplaced into the
      solid, which shows up as dented/wavy flat walls.
- `relax::Union{RelaxOptions, Nothing}`: optional gated vertex-relaxation post-pass
  (Etapa 10), `nothing` = OFF (default). When set, the mesh is relaxed after the final
  connectivity refresh and before export / plane cutting. See [`RelaxOptions`](@ref).
"""
struct MeshGenerationOptions
    scheme::String
    warp_param::Float64
    plane_definitions::Union{Vector{PlaneDefinition},Nothing}
    quality_export::Bool
    cut_points::Symbol
    relax::Union{RelaxOptions,Nothing}

    function MeshGenerationOptions(;
        scheme::String = "A15",
        warp_param::Float64 = 0.3,
        plane_definitions::Union{Vector{PlaneDefinition},Nothing} = nothing,
        quality_export::Bool = false,
        cut_points::Symbol = :linear,
        relax::Union{RelaxOptions,Nothing} = nothing,
    )
        # Validate inputs
        scheme == "A15" || error("Invalid scheme: $scheme. Only 'A15' is supported.")
        warp_param >= 0.0 || error("Invalid warp_param: $warp_param. Must be non-negative.")
        cut_points === :linear || cut_points === :bisection ||
            error("Invalid cut_points: $cut_points. Use :linear or :bisection.")

        new(scheme, warp_param, plane_definitions, quality_export, cut_points, relax)
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

    # Initialize the mesh from structured grid input, then run the shared pipeline.
    mesh = BlockMesh(fine_sdf, fine_grid)
    return generate_tetrahedral_mesh(mesh, output_prefix; options = options)
end

"""
    generate_tetrahedral_mesh(mesh::BlockMesh, output_prefix="output";
                              options=MeshGenerationOptions())

Run the full meshing pipeline on a pre-built [`BlockMesh`](@ref) and return it. This
is the source-agnostic entry point (Etapa 8): it is shared by the structured
file-based method above and is the way to mesh an UNSTRUCTURED HEX8 field --

```julia
src  = simp_to_sdf_source(nodes, hexes, densities)   # or levelset_to_sdf_source(...)
mesh = BlockMesh(src; dx = 1.0)                       # samples a structured lattice
mesh = generate_tetrahedral_mesh(mesh, "part")        # same pipeline, same output contract
```

The pipeline (A15 fill -> warp -> slice -> drop inverted -> keep largest component ->
final connectivity refresh, then optional cutting planes) only reads the field through
`eval_sdf`, so it is identical for both input kinds; only `mesh.sdf_source` differs.
"""
function generate_tetrahedral_mesh(
    mesh::BlockMesh,
    output_prefix::String = "output";
    options::MeshGenerationOptions = MeshGenerationOptions(),
)
    # Generate base tetrahedral mesh using selected discretization scheme
    generate_mesh!(mesh, options.scheme)

    # Warp nodes to isosurface (SDF = 0 level set). The cut-point mode is shared with the
    # slicing below -- the two stages must agree on where the surface crosses an edge.
    warp!(mesh, options.scheme; cut_points = options.cut_points)
    update_connectivity!(mesh; build_ine = false)   # INE not needed until the mesh is final

    # Process isosurface boundary - remove exterior elements
    slice_ambiguous_tetrahedra!(mesh, options.scheme; cut_points = options.cut_points)
    update_connectivity!(mesh; build_ine = false)
    remove_inverted_elements!(mesh)

    # Remove disconnected/isolated mesh components
    remove_isolated_components!(mesh, keep_largest = true)
    update_connectivity!(mesh)                       # final refresh: builds mesh.INE once

    # Optional gated relaxation post-pass (Etapa 10, default OFF). Runs on the final mesh
    # (INE + boundary triangulation valid), before export and before any plane cutting --
    # the plane warp must run last so it is not undone. Topology is unchanged.
    if options.relax !== nothing
        relax_mesh!(mesh, options.relax)
    end

    # Display volume statistics
    @info "Computing mesh volumes..."
    TetMesh_volumes(mesh)

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
