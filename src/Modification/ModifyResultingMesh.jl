# Function to check if a point lies on a bounded plane
function is_on_plane(plane::BoundedPlane, point::Vector{Float64}, tolerance::Float64=1e-10)
  # Check if the point lies in the plane
  dist_to_plane = abs(dot(plane.normal, point - plane.point))
  if dist_to_plane > tolerance
    return false
  end

  # Project the point to the plane's coordinate system
  vec_to_point = point - plane.point

  # Coordinates of the point in the plane
  u_coord = dot(vec_to_point, plane.u)
  v_coord = dot(vec_to_point, plane.v)

  # Check if the projected point lies within the specified shape
  return is_in_shape(plane.shape, u_coord, v_coord)
end

# Check if a point (in plane coordinates) lies within a rectangle
function is_in_shape(shape::Rectangle, u::Float64, v::Float64)
  return abs(u) <= shape.width / 2 && abs(v) <= shape.height / 2
end

# Check if a point (in plane coordinates) lies within a circle
function is_in_shape(shape::Circle, u::Float64, v::Float64)
  return u^2 + v^2 <= shape.radius^2
end

# Check if a point (in plane coordinates) lies within an ellipse
function is_in_shape(shape::Ellipse, u::Float64, v::Float64)
  return (u / shape.a)^2 + (v / shape.b)^2 <= 1
end

function distance_to_bounded_plane(plane::BoundedPlane, point::Vector{Float32})
  # Vector from the plane point to the current point
  vec_to_point = point - plane.point

  # Calculate dot product to determine which side of the plane
  dot_product = dot(plane.normal, vec_to_point)

  # Distance to infinite plane (absolute value for checks)
  dist_to_infinite_plane = abs(dot_product)

  # Project point onto the infinite plane
  projected_point = point - dot_product * plane.normal

  # Check if the projected point lies within the bounded plane
  if is_on_plane(plane, projected_point)
    # If the projection lies within the bounded plane, return distance to the infinite plane
    # with sign - negative in normal direction, positive in opposite direction
    return dot_product > 0 ? -dist_to_infinite_plane : dist_to_infinite_plane
  else
    # If the projection is outside the bounded plane, return a high value
    # while maintaining the correct sign
    return dot_product > 0 ? -1.0e10 : 1.0e10
  end
end

# Function to evaluate the distance of a point from planes (planes_sdf)
function eval_planes_sdf(mesh::BlockMesh, p::SVector{3,Float64}, plane_definitions::Vector{PlaneDefinition})
  # Create bounded planes from definitions
  planes = [BoundedPlane(def.normal, def.point, def.shape) for def in plane_definitions]

  # Initialize with a high positive value
  min_dist = 1.0e10

  # Check distance to each plane
  for plane in planes
    # Convert to Vector{Float32} for compatibility with the distance_to_bounded_plane function
    p_float32 = Vector{Float32}([p[1], p[2], p[3]])
    dist = distance_to_bounded_plane(plane, p_float32)

    # Update minimum distance (comparing absolute values)
    if abs(dist) < abs(min_dist)
      min_dist = dist
    end
  end

  return min_dist
end

# Approximate the gradient of planes_sdf to determine warp direction
function approximate_planes_gradient(mesh::BlockMesh, p::SVector{3,Float64}, plane_definitions::Vector{PlaneDefinition}; h::Float64=1e-3)
  dx = SVector{3,Float64}(h, 0.0, 0.0)
  dy = SVector{3,Float64}(0.0, h, 0.0)
  dz = SVector{3,Float64}(0.0, 0.0, h)

  # Calculate partial derivatives using central differences
  df_dx = (eval_planes_sdf(mesh, p + dx, plane_definitions) - eval_planes_sdf(mesh, p - dx, plane_definitions)) / (2 * h)
  df_dy = (eval_planes_sdf(mesh, p + dy, plane_definitions) - eval_planes_sdf(mesh, p - dy, plane_definitions)) / (2 * h)
  df_dz = (eval_planes_sdf(mesh, p + dz, plane_definitions) - eval_planes_sdf(mesh, p - dz, plane_definitions)) / (2 * h)

  return SVector{3,Float64}(df_dx, df_dy, df_dz)
end


# Function to move a node to the zero level of planes_sdf
function warp_node_to_planes_isocontour!(mesh::BlockMesh, node_index::Int, plane_definitions::Vector{PlaneDefinition}, max_iter::Int)
  tol = mesh.grid_tol
  current_position = mesh.X[node_index]

  for iter in 1:max_iter
    # Evaluate planes_sdf at current position
    f = eval_planes_sdf(mesh, current_position, plane_definitions)

    # If we're close enough to the isosurface, end
    abs2(f) < tol * tol && break

    # Calculate gradient for displacement direction
    grad = approximate_planes_gradient(mesh, current_position, plane_definitions)
    norm_grad_squared = sum(abs2, grad)

    # If gradient is too small, end
    norm_grad_squared < 1e-16 && break

    # Newton step
    dp = (f / norm_grad_squared) * grad
    current_position -= dp
  end

  # Calculate current planes_sdf value after warping
  current_sdf = eval_planes_sdf(mesh, current_position, plane_definitions)

  # Update node position and its SDF value
  mesh.X[node_index] = current_position
  if abs(current_sdf) < tol * 2
    mesh.node_sdf[node_index] = 0.0
  else
    mesh.node_sdf[node_index] = current_sdf
  end
end


function warp_node_to_planes_isocontour_new!(mesh::BlockMesh, node_index::Int, plane_definitions::Vector{PlaneDefinition}, max_iter::Int)
  # Use a tighter tolerance for convergence
  tol = mesh.grid_tol * 0.01  # Much stricter tolerance
  current_position = mesh.X[node_index]

  # Track the best position and its SDF value
  best_position = copy(current_position)
  best_sdf = eval_planes_sdf(mesh, current_position, plane_definitions)
  best_sdf_abs = abs(best_sdf)

  # Use adaptive step size with line search
  for iter in 1:max_iter
    # Evaluate planes_sdf at current position
    f = eval_planes_sdf(mesh, current_position, plane_definitions)

    # If we're close enough to the isosurface, end
    if abs(f) < tol
      best_position = current_position
      best_sdf = f
      break
    end

    # Calculate gradient for displacement direction
    grad = approximate_planes_gradient(mesh, current_position, plane_definitions)
    grad_norm = sqrt(sum(abs2, grad))

    # If gradient is too small, end
    if grad_norm < 1e-10
      break
    end

    # Normalize gradient
    normalized_grad = grad / grad_norm

    # Base step size - adaptive based on distance
    step_size = min(abs(f), 0.1 * mesh.grid_tol)

    # Line search with multiple step sizes
    for alpha in [1.0, 0.5, 0.25, 0.1, 0.01]
      # Calculate new position with current step size
      test_position = current_position - alpha * step_size * normalized_grad

      # Evaluate SDF at test position
      test_sdf = eval_planes_sdf(mesh, test_position, plane_definitions)
      test_sdf_abs = abs(test_sdf)

      # If this position is better (closer to zero), use it
      if test_sdf_abs < best_sdf_abs
        best_position = test_position
        best_sdf = test_sdf
        best_sdf_abs = test_sdf_abs

        # If we're very close to the surface, stop line search
        if test_sdf_abs < tol
          break
        end
      end
    end

    # Update current position to best found
    current_position = best_position

    # If we're close enough, stop iterations
    if best_sdf_abs < tol
      break
    end

    # If we aren't making progress anymore, stop
    if iter > 5 && best_sdf_abs > 0.99 * abs(eval_planes_sdf(mesh, mesh.X[node_index], plane_definitions))
      break
    end
  end

  # Use direct projection for final precision if we have a plane surface
  # This is especially useful for planar surfaces where we can project directly
  for def in plane_definitions
    # For efficiency, only check planes that are close enough
    if abs(dot(def.normal, best_position - def.point)) < 10 * tol
      # Project the point directly onto the plane
      projected_point = best_position - dot(def.normal, best_position - def.point) * def.normal

      # Check if the projected point lies within the bounded plane
      if is_on_plane(BoundedPlane(def.normal, def.point, def.shape), projected_point)
        # Use this exact projection
        best_position = projected_point
        best_sdf = 0.0  # Exactly on the plane
        break
      end
    end
  end

  # Update node position and its SDF value
  mesh.X[node_index] = best_position
  mesh.node_sdf[node_index] = best_sdf

  # Ensure the SDF is exactly zero if we're very close to the plane
  if abs(best_sdf) < tol
    mesh.node_sdf[node_index] = 0.0
  end

  # No return value - function modifies mesh in-place
end


# Main function for modifying the mesh according to planes_sdf
function warp_mesh_by_planes_sdf!(mesh::BlockMesh, plane_definitions::Vector{PlaneDefinition}, warp_param::Float64; max_iter::Int=20)
  # Check if any plane definitions were provided
  if isempty(plane_definitions)
    @info "No plane definitions provided, skipping plane-based warping."
    return
  end

  # Calculate the longest edge and then the threshold for displacement - same logic as in the warp! function
  threshold_sdf = warp_param * mesh.grid_step

  @info "Planes warping: max edge = $max_edge, threshold_sdf = $threshold_sdf"

  # Preparation for warping - find nodes near planes (in one direction)
  nodes_to_warp = Int[]
  for i in 1:length(mesh.X)
    plane_sdf = eval_planes_sdf(mesh, mesh.X[i], plane_definitions)
    # If the node is close enough to the plane (from both sides), add it to nodes to warp
    # if abs(plane_sdf) < threshold_sdf
    # if plane_sdf < threshold_sdf
    if plane_sdf < 0. # one-side warp
      push!(nodes_to_warp, i)
    end
  end

  @info "Found $(length(nodes_to_warp)) nodes near planes for warping"

  # Warp nodes near planes to the zero level
  for node_idx in nodes_to_warp
    warp_node_to_planes_isocontour!(mesh, node_idx, plane_definitions, max_iter)
  end

  # Update mesh topology
  update_connectivity!(mesh)

  # Remove elements that have all nodes at the zero level or don't have at least one node with positive SDF
  new_IEN = Vector{Vector{Int64}}()
  for tet in mesh.IEN
    # Get SDF values for all nodes of the tetrahedron
    tet_sdf = [eval_planes_sdf(mesh, mesh.X[i], plane_definitions) for i in tet]

    # Count nodes at zero level
    zero_nodes = count(x -> abs(x) < mesh.grid_tol, tet_sdf)

    # Count nodes with positive SDF value
    pos_nodes = count(x -> x > 0, tet_sdf)

    # Keep the tetrahedron only if:
    # 1. It has at least one node with positive value (relative to planes_sdf)
    # 2. It doesn't have all nodes at the zero level
    if pos_nodes > 0 && zero_nodes < 4
      push!(new_IEN, tet)
    end
  end

  # Update connectivity
  mesh.IEN = new_IEN
  @info "After removing unsuitable elements: $(length(mesh.IEN)) tetrahedra"

  # Remove nodes outside the isocontour and elements that contain them
  # remove_nodes_outside_isocontour!(mesh)

  # Remove any inverted elements (with negative Jacobian determinant)
  remove_inverted_elements!(mesh)

  # Final topology update
  update_connectivity!(mesh)

  @info "Mesh modification according to planes_sdf completed"
end

#TODO: Create a mesh trimming that will support element cutting and subsequent mesh optimization.

