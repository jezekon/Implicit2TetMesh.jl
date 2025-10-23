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

    function BlockMesh(fine_sdf::Array, fine_grid::Array)

        # Convert fine_grid into a 3D array of SVectors
        grid = Array{SVector{3,Float64},3}(undef, size(fine_grid))
        for i in eachindex(fine_grid)
            # Assume fine_grid[i] is an array of Float64; convert to SVector
            grid[i] = SVector{3,Float64}(fine_grid[i]...)
        end
        step = maximum(abs.(grid[1, 1, 1] - grid[2, 2, 2]))
        sdf = Float64.(fine_sdf)
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

        return mesh
    end
end

