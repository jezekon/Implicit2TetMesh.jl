# ----------------------------
# Helper function: Get SDF values at the 8 corners of a cell (unchanged)
# ----------------------------
function get_cell_sdf_values(mesh::BlockMesh, i::Int, j::Int, k::Int)
  1 <= i < mesh.nx || throw(BoundsError(mesh.SDF, i))
  1 <= j < mesh.ny || throw(BoundsError(mesh.SDF, j))
  1 <= k < mesh.nz || throw(BoundsError(mesh.SDF, k))
  return SVector{8}(
    mesh.SDF[i, j, k],       # front-bottom-left
    mesh.SDF[i+1, j, k],     # front-bottom-right
    mesh.SDF[i+1, j+1, k],   # front-top-right
    mesh.SDF[i, j+1, k],     # front-top-left
    mesh.SDF[i, j, k+1],     # back-bottom-left
    mesh.SDF[i+1, j, k+1],   # back-bottom-right
    mesh.SDF[i+1, j+1, k+1], # back-top-right
    mesh.SDF[i, j+1, k+1]    # back-top-left
  )
end


# ----------------------------
# Helper function: SDF evaluation using trilinear interpolation
# ----------------------------
function eval_sdf(mesh::BlockMesh, p::SVector{3,Float64})
  # Get minimum and maximum grid coordinates
  vmin = mesh.grid[1, 1, 1]
  vmax = mesh.grid[end, end, end]
  
  # Normalize point coordinates to interval [0,1]
  r = (p .- vmin) ./ (vmax .- vmin)
  
  # Convert to grid indices
  i_f = r[1] * (mesh.nx - 1) + 1
  j_f = r[2] * (mesh.ny - 1) + 1
  k_f = r[3] * (mesh.nz - 1) + 1
  
  i0 = clamp(floor(Int, i_f), 1, mesh.nx - 1)
  j0 = clamp(floor(Int, j_f), 1, mesh.ny - 1)
  k0 = clamp(floor(Int, k_f), 1, mesh.nz - 1)
  
  i1 = i0 + 1
  j1 = j0 + 1
  k1 = k0 + 1
  
  # Local weights
  xd = i_f - i0
  yd = j_f - j0
  zd = k_f - k0
  
  # Get SDF values at the eight corners of the cell
  c000 = mesh.SDF[i0, j0, k0]
  c100 = mesh.SDF[i1, j0, k0]
  c010 = mesh.SDF[i0, j1, k0]
  c110 = mesh.SDF[i1, j1, k0]
  c001 = mesh.SDF[i0, j0, k1]
  c101 = mesh.SDF[i1, j0, k1]
  c011 = mesh.SDF[i0, j1, k1]
  c111 = mesh.SDF[i1, j1, k1]
  
  # Trilinear interpolation
  c00 = c000 * (1 - xd) + c100 * xd
  c01 = c001 * (1 - xd) + c101 * xd
  c10 = c010 * (1 - xd) + c110 * xd
  c11 = c011 * (1 - xd) + c111 * xd
  
  c0 = c00 * (1 - yd) + c10 * yd
  c1 = c01 * (1 - yd) + c11 * yd
  
  return c0 * (1 - zd) + c1 * zd
end


# ----------------------------
# Helper function: Approximation of the SDF gradient at point p using central differences
# ----------------------------
function approximate_gradient(mesh::BlockMesh, p::SVector{3,Float64}; h::Float64=1e-3)
  dx = SVector{3,Float64}(h, 0.0, 0.0)
  dy = SVector{3,Float64}(0.0, h, 0.0)
  dz = SVector{3,Float64}(0.0, 0.0, h)
  df_dx = (eval_sdf(mesh, p + dx) - eval_sdf(mesh, p - dx)) / (2 * h)
  df_dy = (eval_sdf(mesh, p + dy) - eval_sdf(mesh, p - dy)) / (2 * h)
  df_dz = (eval_sdf(mesh, p + dz) - eval_sdf(mesh, p - dz)) / (2 * h)
  return SVector{3,Float64}(df_dx, df_dy, df_dz)
end


# Estimate the gradient of the SDF around the node
# First, let's modify compute_gradient to work with both node indices and positions
function compute_gradient(mesh::BlockMesh, p::SVector{3,Float64}; δ::Float64=1e-3)
  # Pre-allocate unit vectors as static vectors for better performance
  unit_vectors = (
    SVector{3,Float64}(1.0, 0.0, 0.0),
    SVector{3,Float64}(0.0, 1.0, 0.0),
    SVector{3,Float64}(0.0, 0.0, 1.0)
  )

  # Use tuple comprehension for better compile-time optimization
  grad = ntuple(3) do d
    unit_vec = unit_vectors[d]
    # Compute central difference
    (eval_sdf(mesh, p + δ * unit_vec) - eval_sdf(mesh, p - δ * unit_vec)) / (2δ)
  end

  return SVector{3,Float64}(grad)
end

# Add method for node index input
function compute_gradient(mesh::BlockMesh, node_index::Int; δ::Float64=1e-3)
  return compute_gradient(mesh, mesh.X[node_index]; δ=δ)
end

