# ----------------------------
# Quantization function – unchanged
# ----------------------------
function quantize(p::SVector{3,Float64}, tol::Float64)
  return (round(p[1] / tol) * tol, round(p[2] / tol) * tol, round(p[3] / tol) * tol)
end

# ----------------------------
# Function for discretizing a cell using A15 scheme
# ----------------------------
function process_cell_A15!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  tol = mesh.grid_tol
  
  # First check the SDF values of the current cell
  current_sdf_values = get_cell_sdf_values(mesh, i, j, k)
  
  # If any value in the current cell is positive or close to zero, we definitely process the cell
  if any(x -> x >= -tol, current_sdf_values)
    # Continue with regular processing
  else
    # All values in the current cell are negative - check neighboring cells
  
    # Define offsets for neighboring cells (direct neighbors in all directions)
    neighbor_offsets = [
      (1,0,0), (-1,0,0),   # neighbors in x direction
      (0,1,0), (0,-1,0),   # neighbors in y direction
      (0,0,1), (0,0,-1)    # neighbors in z direction
    ] #TODO: select only relevant one
    
    # Initialize a flag indicating that all neighboring cells have negative SDF values
    all_neighbors_negative = true
    
    # Check SDF values of neighboring cells
    for (di, dj, dk) in neighbor_offsets
      # Calculate indices of neighboring cell
      ni, nj, nk = i + di, j + dj, k + dk
      
      # Check grid boundaries
      if 1 <= ni < mesh.nx && 1 <= nj < mesh.ny && 1 <= nk < mesh.nz
        # Get SDF values for neighboring cell
        neighbor_sdf = get_cell_sdf_values(mesh, ni, nj, nk)
        
        # If any value in the neighboring cell has a positive or zero value,
        # set the flag and end the check
        if any(x -> x >= -tol, neighbor_sdf)
          all_neighbors_negative = false
          break
        end
      end
    end
        
    # If all values in the current cell and all neighboring cells are negative,
    # we can safely skip the cell
    if all_neighbors_negative
      return
    end
  end

  # Retrieve min and max corners of the cell
  v000 = mesh.grid[i, j, k]
  v111 = mesh.grid[i+1, j+1, k+1]
  vmins = v000
  vmaxs = v111

  # Precompute differences for coordinate interpolation
  Δ = vmaxs .- vmins
  local_mapping = Dict{Int,Int}()

  @inbounds for li in 1:length(tile_ref)
    # tile_ref is assumed to be defined in A15_scheme.jl and normalized (in [0,1]*4 originally)
    local_coord = SVector{3,Float64}(tile_ref[li] ./ 4.0)  # normalized coordinates in [0,1]
    # Compute physical point using linear interpolation
    p = vmins .+ Δ .* local_coord
    p = SVector{3,Float64}(p)  # ensure static vector type
    p_key = quantize(p, tol)
    if haskey(mesh.node_hash, p_key)
      local_mapping[li] = mesh.node_hash[p_key]
    else
      push!(mesh.X, p)

      sdf_aprox = eval_sdf(mesh, p)
      push!(mesh.node_sdf, sdf_aprox)
      local_index = length(mesh.X)
      local_mapping[li] = local_index
      mesh.node_hash[p_key] = local_index
    end
  end

  # Process tetrahedral connectivity from A15 scheme
  @inbounds for tet in tetra_connectivity
    global_tet = [local_mapping[li] for li in tet]
  
    # Get the coordinates of the tetrahedron vertices
    tet_coords = [mesh.X[idx] for idx in global_tet]
  
    # Directly evaluate SDF at each vertex position for maximum accuracy
    # This is more accurate than using pre-computed values
    tet_sdf = [eval_sdf(mesh, coord) for coord in tet_coords]
  
    # Update the stored SDF values with these more accurate evaluations
    for (i, idx) in enumerate(global_tet)
      mesh.node_sdf[idx] = tet_sdf[i]
    end
  
    # Include tetrahedron only if at least one vertex is inside or on the boundary
    if any(x -> x >= 0, tet_sdf)
      push!(mesh.IEN, global_tet)
    end
  end
end

# ----------------------------
# Function for discretizing a cell using Schlafli orthoscheme (unchanged logic, only minor type annotation changes)
# ----------------------------
function process_cell_Schlafli!(mesh::BlockMesh, i::Int, j::Int, k::Int)
  tol = mesh.grid_tol
  # Get SDF values at the 8 corners of the cell
  sdf_values = get_cell_sdf_values(mesh, i, j, k)
  if !any(x -> x >= -tol, sdf_values)
    return
  end

  local_mapping = Dict{Int,Int}()
  # Define cell nodes as SVectors from grid
  cell_nodes = [
    mesh.grid[i, j, k],     # Node 1: front-bottom-left
    mesh.grid[i+1, j, k],       # Node 2: front-bottom-right
    mesh.grid[i+1, j+1, k],       # Node 3: front-top-right
    mesh.grid[i, j+1, k],       # Node 4: front-top-left
    mesh.grid[i, j, k+1],     # Node 5: back-bottom-left
    mesh.grid[i+1, j, k+1],     # Node 6: back-bottom-right
    mesh.grid[i+1, j+1, k+1],     # Node 7: back-top-right
    mesh.grid[i, j+1, k+1]      # Node 8: back-top-left
  ]

  @inbounds for li in 1:8
    p = cell_nodes[li]
    p_key = quantize(p, tol)
    if haskey(mesh.node_hash, p_key)
      local_mapping[li] = mesh.node_hash[p_key]
    else
      push!(mesh.X, p)
      push!(mesh.node_sdf, sdf_values[li])
      local_index = length(mesh.X)
      local_mapping[li] = local_index
      mesh.node_hash[p_key] = local_index
    end
  end

  # Construct tetrahedra according to Schlafli scheme
  @inbounds for tet in schlafli_tet_connectivity
    global_tet = [local_mapping[li] for li in tet]
    tet_sdf = [mesh.node_sdf[idx] for idx in global_tet]
    if any(x -> x >= 0, tet_sdf)
      push!(mesh.IEN, global_tet)
    end
  end
end

# ----------------------------
# Merge duplicate nodes after mesh generation (unchanged logic, only type annotations updated)
# ----------------------------
function merge_duplicate_nodes!(mesh::BlockMesh)
  @info "Merging duplicate nodes"
  println("  Before merging duplicates: $(length(mesh.X)) nodes")
  tol = mesh.grid_tol*10^4
  new_nodes = Vector{SVector{3,Float64}}()
  new_node_sdf = Vector{Float64}()
  node_map = Dict{Int,Int}()
  global_hash = Dict{NTuple{3,Float64},Int}()
  @inbounds for i in 1:length(mesh.X)
    p = mesh.X[i]
    p_key = quantize(p, tol)
    if haskey(global_hash, p_key)
      node_map[i] = global_hash[p_key]
    else
      push!(new_nodes, p)
      push!(new_node_sdf, mesh.node_sdf[i])
      new_index = length(new_nodes)
      node_map[i] = new_index
      global_hash[p_key] = new_index
    end
  end
  @inbounds for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [node_map[old] for old in mesh.IEN[i]]
  end
  mesh.X = new_nodes
  mesh.node_sdf = new_node_sdf
  mesh.node_map = node_map
  println("  After merging duplicates:  $(length(mesh.X)) nodes")
end

# ----------------------------
# Cleanup unused nodes and reindex connectivity (unchanged logic)
# ----------------------------
function cleanup_unused_nodes!(mesh::BlockMesh)
  @info "Cleaning unused nodes"
  println("  Number of nodes before cleanup: $(length(mesh.X))")
  used_nodes = Set{Int64}()
  @inbounds for element in mesh.IEN
    union!(used_nodes, element)
  end
  new_node_map = Dict{Int64,Int64}()
  new_coords = Vector{SVector{3,Float64}}()
  new_node_sdf = Vector{Float64}()
  sorted_used = sort(collect(used_nodes))
  for (new_id, old_id) in enumerate(sorted_used)
    new_node_map[old_id] = new_id
    push!(new_coords, mesh.X[old_id])
    push!(new_node_sdf, mesh.node_sdf[old_id])
  end
  @inbounds for i in 1:length(mesh.IEN)
    mesh.IEN[i] = [new_node_map[old_id] for old_id in mesh.IEN[i]]
  end
  mesh.X = new_coords
  mesh.node_sdf = new_node_sdf
  mesh.node_map = new_node_map
  println("  Number of nodes after cleanup:  $(length(mesh.X))")
end

# ----------------------------
# Create inverse connectivity (unchanged logic)
# ----------------------------
function create_INE!(mesh::BlockMesh)
  mesh.INE = [Vector{Int64}() for _ in 1:length(mesh.X)]
  @inbounds for (elem_id, element) in enumerate(mesh.IEN)
    for node_id in element
      push!(mesh.INE[node_id], elem_id)
    end
  end
  return mesh
end

# ----------------------------
# Modified mesh generation function with scheme selection
# ----------------------------
function generate_mesh!(mesh::BlockMesh, scheme::String)
  empty!(mesh.X)
  empty!(mesh.IEN)
  empty!(mesh.node_sdf)
  empty!(mesh.node_map)
  empty!(mesh.cell_center_map)
  empty!(mesh.node_hash)

  @info "Generating mesh with $(scheme) scheme..."

  @inbounds for i in 1:mesh.nx-1
    for j in 1:mesh.ny-1
      for k in 1:mesh.nz-1
        if scheme == "A15"
          process_cell_A15!(mesh, i, j, k)
        elseif scheme == "Schlafli"
          process_cell_Schlafli!(mesh, i, j, k)
        else
          error("Unknown scheme: $scheme")
        end
      end
    end
  end

  cleanup_unused_nodes!(mesh)
  create_INE!(mesh)
end


# Function to calculate the length of the longest edge between nodes in all tetrahedra
# TODO: Just take the initial discretization and take only one element (for speedup)
function longest_edge(mesh::BlockMesh)
  # Pre-allocate maximum length with type stability
  max_length = zero(eltype(mesh.X[1]))

  # Vectorize operations by using broadcast
  for tet in mesh.IEN
    # Use array views for better memory efficiency
    nodes = @view mesh.X[tet]
    # Calculate all edge lengths at once using comprehension
    lengths = [norm(nodes[i] - nodes[j]) for i in 1:3 for j in i+1:4]
    # Use built-in maximum function
    max_length = max(max_length, maximum(lengths))
  end
  return max_length
end


# Update warp_node_to_isocontour! to handle positions directly
function warp_node_to_isocontour!(mesh::BlockMesh, node_index::Int, max_dist::Float64, max_iter)
  tol = mesh.grid_tol
  current_position = mesh.X[node_index]

  for iter in 1:max_iter
    f = eval_sdf(mesh, current_position)

    # Early return if we're close enough to the isocontour
    abs2(f) < tol * tol && break

    grad = compute_gradient(mesh, current_position)
    norm_grad_squared = sum(abs2, grad)

    # Early return if gradient is too small
    norm_grad_squared < 1e-16 && break

    # Newton step
    dp = (f / norm_grad_squared) * grad
    current_position -= dp
  end

  norm_dist = norm(current_position .- mesh.X[node_index])

  if norm_dist <= (max_dist * 2)
    current_sdf = eval_sdf(mesh, current_position)
    if abs(current_sdf) < tol*4
      mesh.node_sdf[node_index] = 0.
    else
      # println("current_sdf: ", current_sdf)
      mesh.node_sdf[node_index] = current_sdf
    end
    mesh.X[node_index] = current_position
  end
end

# Main function for node warping - ordered warping
#
# First, nodes with positive SDF value (inside the isosurface) are adjusted, 
# then nodes with negative values.
# Nodes are moved toward the zero level of SDF (isosurface) and the displacement threshold 
# is calculated as threshold_sdf = 0.5 * (length of the longest tetrahedral edge).
function warp!(mesh::BlockMesh, scheme::String, max_iter::Int=160)
  # Calculate the longest edge and then the threshold for displacement
  @info "Warping nodes to isocontour..."
  max_edge = longest_edge(mesh)
  if scheme == "A15"
    threshold_sdf = 0.15 * mesh.grid_step
  elseif scheme =="Schlafli"
    threshold_sdf = 0.3 * mesh.grid_step
  else
    @error "Unknown scheme"
  end

  println("  max edge length = $(round(max_edge, sigdigits=4)) , threshold sdf for warp = $threshold_sdf")
  # First pass: nodes with positive SDF value (inside)
  for i in 1:length(mesh.X)
    sdf = mesh.node_sdf[i]
    if sdf > 0 && abs(sdf) < threshold_sdf
      warp_node_to_isocontour!(mesh, i, threshold_sdf, max_iter)
    end
  end
  # Second pass: nodes with negative SDF value (outside)
  for i in 1:length(mesh.X)
    sdf = mesh.node_sdf[i]
    if sdf < 0 && abs(sdf) < threshold_sdf
      warp_node_to_isocontour!(mesh, i, threshold_sdf, max_iter)
    end
  end
  #TODO:Warp only one node from element wich is close to boundary (now it warps every node which is close enough)
end

# ---------------------------------------------------
# Function: Update mesh topology (mesh.X, mesh.IEN, mesh.INE)
# ---------------------------------------------------
function update_connectivity!(mesh::BlockMesh)
  cleanup_unused_nodes!(mesh)        # Recalculates mesh.X, mesh.node_sdf and reindexes mesh.IEN and mesh.node_map
  merge_duplicate_nodes!(mesh)       # Merges duplicate nodes and adjusts connectivity in mesh.IEN
  create_INE!(mesh)                  # Creates inverse connectivity (mesh.INE)
end
