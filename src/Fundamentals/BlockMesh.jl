# ----------------------------
# Block mesh structure
# ----------------------------
mutable struct BlockMesh
    nx::Int                                    # Number of grid points in x-direction
    ny::Int                                    # Number of grid points in y-direction
    nz::Int                                    # Number of grid points in z-direction
    grid::Array{SVector{3,Float64},3}          # 3D array of node coordinates (basic grid) using static vectors
    grid_step::Float64                         # Spatial step size between adjacent grid points (uniform spacing)
    grid_tol::Float64                          # Geometric tolerance
    SDF::Array{Float64,3}                      # Signed Distance Function values at each grid point
    X::Vector{SVector{3,Float64}}              # List of physical node coordinates (nodes used in mesh)
    IEN::Vector{Vector{Int64}}                 # Tetrahedral connectivity (elements)
    INE::Vector{Vector{Int64}}                 # Inverse connectivity: for each node, list of adjacent elements
    node_sdf::Vector{Float64}                  # SDF values at nodes
    node_map::Dict{Int64,Int64}                # Mapping from original grid node index -> new node index
    cell_center_map::Dict{Tuple{Int,Int,Int},Int64}  # Mapping for cell centers (Steiner points)
    node_hash::Dict{NTuple{3,Float64},Int64}   # Global dictionary for merging nodes
    sdf_source::SDFSource                       # Pluggable field source; mesh.SDF is its lattice-corner cache

    function BlockMesh(fine_sdf::Array, fine_grid::Array)

        # Convert fine_grid into a 3D array of SVectors
        grid = Array{SVector{3,Float64},3}(undef, size(fine_grid))
        for i in eachindex(fine_grid)
            # Assume fine_grid[i] is an array of Float64; convert to SVector
            grid[i] = SVector{3,Float64}(fine_grid[i]...)
        end
        step = maximum(abs.(grid[1, 1, 1] - grid[2, 2, 2]))
        # Convention: phi < 0 = INSIDE, phi > 0 = OUTSIDE, phi = 0 = ON the surface,
        # matching the reference implementations (quartet / Labelle / isostuffer).
        # The input data files use the opposite sign (positive = inside), so we negate
        # the field once here. This is the single source for mesh.SDF; everything else
        # reads it through get_cell_sdf_values / eval_sdf. The data files on disk are
        # left untouched.
        sdf = -Float64.(fine_sdf)
        nx, ny, nz = size(grid)

        mesh = new(nx, ny, nz)
        mesh.grid = grid
        mesh.grid_step = step
        mesh.grid_tol = 1e-8 * step
        mesh.SDF = sdf
        mesh.node_map = Dict{Int64,Int64}()
        mesh.cell_center_map = Dict{Tuple{Int,Int,Int},Int64}()
        mesh.X = Vector{SVector{3,Float64}}()
        mesh.IEN = Vector{Vector{Int64}}()
        mesh.INE = Vector{Vector{Int64}}()
        mesh.node_sdf = Vector{Float64}()
        mesh.node_hash = Dict{NTuple{3,Float64},Int64}()
        # Structured input: the generation lattice IS the input grid, and the source
        # wraps the SAME grid + negated values. eval_sdf therefore returns exactly the
        # original trilinear interpolation (bit-identical), while get_cell_sdf_values
        # keeps reading mesh.SDF (= the source's values array).
        mesh.sdf_source = StructuredSDF(grid, sdf)

        return mesh
    end

    # Generic constructor: build the structured generation lattice from an arbitrary
    # SDF source (Etapa 8). Used for unstructured HEX8 inputs, whose source has no
    # native grid. See the keyword constructor below for documentation.
    function BlockMesh(source::SDFSource, grid::Array{SVector{3,Float64},3},
                       sdf::Array{Float64,3}, step::Float64)
        nx, ny, nz = size(grid)
        mesh = new(nx, ny, nz)
        mesh.grid = grid
        mesh.grid_step = step
        mesh.grid_tol = 1e-8 * step
        mesh.SDF = sdf
        mesh.node_map = Dict{Int64,Int64}()
        mesh.cell_center_map = Dict{Tuple{Int,Int,Int},Int64}()
        mesh.X = Vector{SVector{3,Float64}}()
        mesh.IEN = Vector{Vector{Int64}}()
        mesh.INE = Vector{Vector{Int64}}()
        mesh.node_sdf = Vector{Float64}()
        mesh.node_hash = Dict{NTuple{3,Float64},Int64}()
        mesh.sdf_source = source
        return mesh
    end
end

"""
    BlockMesh(source::SDFSource; dx, padding = 2)

Build a `BlockMesh` whose structured generation lattice is sampled from an arbitrary
[`SDFSource`](@ref). This is the entry point for unstructured HEX8 inputs (Etapa 8):
the mesher still fills a STRUCTURED lattice (a property of isosurface stuffing), but
the field at the lattice corners -- and at every later `eval_sdf` query -- comes from
the source.

The lattice spans `bbox(source)` grown by a `padding`-cell ring on every side
(mirroring quartet's `main.cpp` grid sizing), with uniform spacing `dx`. Each lattice
corner caches `eval_sdf(source, corner)` into `mesh.SDF`; the padding ring sits
outside the source's domain and so reads as exterior (positive), keeping the zero
isosurface fully enclosed.

# Arguments
- `dx::Real`: lattice spacing (cell edge length). Required.
- `padding::Integer = 2`: number of extra lattice cells added around the bounding box
  on each side.
"""
function BlockMesh(source::SDFSource; dx::Real, padding::Integer = 2)
    dx > 0 || error("BlockMesh: dx must be positive, got $dx")
    padding >= 0 || error("BlockMesh: padding must be non-negative, got $padding")
    dxf = Float64(dx)

    (bmin, bmax) = bbox(source)
    origin = bmin .- padding * dxf

    # Number of corner points so the lattice covers bmax + padding*dx on each axis.
    span = (bmax .+ padding * dxf) .- origin
    nx = ceil(Int, span[1] / dxf) + 1
    ny = ceil(Int, span[2] / dxf) + 1
    nz = ceil(Int, span[3] / dxf) + 1

    # Build the lattice node coordinates and sample the source at each corner.
    grid = Array{SVector{3,Float64},3}(undef, nx, ny, nz)
    sdf = Array{Float64,3}(undef, nx, ny, nz)
    for k = 1:nz
        for j = 1:ny
            for i = 1:nx
                p = SVector{3,Float64}(
                    origin[1] + (i - 1) * dxf,
                    origin[2] + (j - 1) * dxf,
                    origin[3] + (k - 1) * dxf,
                )
                grid[i, j, k] = p
                sdf[i, j, k] = eval_sdf(source, p)
            end
        end
    end

    return BlockMesh(source, grid, sdf, dxf)
end

