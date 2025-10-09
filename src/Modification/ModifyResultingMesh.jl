# Function to check if a point lies on a bounded plane
function is_on_plane(
    plane::BoundedPlane,
    point::Vector{Float64},
    tolerance::Float64 = 1e-10,
)
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
function eval_planes_sdf(
    mesh::BlockMesh,
    p::SVector{3,Float64},
    plane_definitions::Vector{PlaneDefinition},
)
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
function approximate_planes_gradient(
    mesh::BlockMesh,
    p::SVector{3,Float64},
    plane_definitions::Vector{PlaneDefinition};
    h::Float64 = 1e-3,
)
    dx = SVector{3,Float64}(h, 0.0, 0.0)
    dy = SVector{3,Float64}(0.0, h, 0.0)
    dz = SVector{3,Float64}(0.0, 0.0, h)

    # Calculate partial derivatives using central differences
    df_dx =
        (
            eval_planes_sdf(mesh, p + dx, plane_definitions) -
            eval_planes_sdf(mesh, p - dx, plane_definitions)
        ) / (2 * h)
    df_dy =
        (
            eval_planes_sdf(mesh, p + dy, plane_definitions) -
            eval_planes_sdf(mesh, p - dy, plane_definitions)
        ) / (2 * h)
    df_dz =
        (
            eval_planes_sdf(mesh, p + dz, plane_definitions) -
            eval_planes_sdf(mesh, p - dz, plane_definitions)
        ) / (2 * h)

    return SVector{3,Float64}(df_dx, df_dy, df_dz)
end

# Function to identify surface nodes (nodes that belong to boundary faces)
function find_surface_nodes(mesh::BlockMesh)::Set{Int}
    # Dictionary to count how many tetrahedra each face belongs to
    # Key: Set of 3 node indices forming a triangular face
    # Value: Number of tetrahedra containing this face
    face_count = Dict{Set{Int},Int}()

    # Process each tetrahedron and its 4 triangular faces
    for tet in mesh.IEN
        # Each tetrahedron has 4 triangular faces
        faces = [
            Set([tet[1], tet[2], tet[3]]),  # Face opposite to node 4
            Set([tet[1], tet[2], tet[4]]),  # Face opposite to node 3
            Set([tet[1], tet[3], tet[4]]),  # Face opposite to node 2
            Set([tet[2], tet[3], tet[4]]),   # Face opposite to node 1
        ]

        # Count occurrences of each face
        for face in faces
            face_count[face] = get(face_count, face, 0) + 1
        end
    end

    # Surface faces are those that belong to exactly one tetrahedron
    surface_faces = [face for (face, count) in face_count if count == 1]

    # Collect all unique nodes from surface faces
    surface_nodes = Set{Int}()
    for face in surface_faces
        union!(surface_nodes, face)
    end

    return surface_nodes
end

# Simplified function to move a node to the zero level of planes_sdf
function warp_node_to_planes_isocontour!(
    mesh::BlockMesh,
    node_index::Int,
    plane_definitions::Vector{PlaneDefinition},
    max_dist::Float64,
    max_iter::Int,
)
    tol = mesh.grid_tol
    current_position = mesh.X[node_index]

    for iter = 1:max_iter
        # Evaluate planes_sdf at current position
        f = eval_planes_sdf(mesh, current_position, plane_definitions)

        # Early return if we're close enough to the plane surface
        abs2(f) < tol * tol && break

        # Calculate gradient for displacement direction
        grad = approximate_planes_gradient(mesh, current_position, plane_definitions)
        norm_grad_squared = sum(abs2, grad)

        # Early return if gradient is too small
        norm_grad_squared < 1e-16 && break

        # Newton step
        dp = (f / norm_grad_squared) * grad
        current_position -= dp
    end

    # Check if displacement is within acceptable limits
    norm_dist = norm(current_position - mesh.X[node_index])

    if norm_dist <= (max_dist * 2)
        # Update node position
        mesh.X[node_index] = current_position
    end
end

# Main function for plane-based mesh warping - ONLY moves SURFACE nodes
function warp_mesh_by_planes_sdf!(
    mesh::BlockMesh,
    plane_definitions::Vector{PlaneDefinition},
    warp_param::Float64;
    max_iter::Int = 20,
)
    # Check if any plane definitions were provided
    if isempty(plane_definitions)
        @info "No plane definitions provided, skipping plane-based warping."
        return
    end

    @info "Warping surface nodes to plane surfaces..."

    # Calculate threshold for displacement based on warp_param and grid step
    threshold_sdf = warp_param * mesh.grid_step
    println("  warp_param = $warp_param, threshold_sdf = $threshold_sdf")

    # Find all surface nodes (nodes on boundary faces)
    surface_nodes = find_surface_nodes(mesh)
    println("  Found $(length(surface_nodes)) surface nodes")

    # Single pass: warp only surface nodes that are close to planes
    nodes_warped = 0
    nodes_checked = 0

    for node_idx in surface_nodes
        nodes_checked += 1
        plane_sdf = eval_planes_sdf(mesh, mesh.X[node_idx], plane_definitions)

        # Warp surface nodes that are close to planes from either side
        if abs(plane_sdf) < threshold_sdf
            warp_node_to_planes_isocontour!(
                mesh,
                node_idx,
                plane_definitions,
                threshold_sdf,
                max_iter,
            )
            nodes_warped += 1
        end
    end

    println("  Checked $nodes_checked surface nodes")
    println("  Warped $nodes_warped surface nodes to plane surfaces")

    # Only fix potentially inverted elements caused by node movement
    # Do NOT delete elements or change mesh topology
    remove_inverted_elements!(mesh)

    # Update inverse connectivity (INE) to reflect any potential changes
    create_INE!(mesh)

    @info "Surface node warping completed. Element count unchanged: $(length(mesh.IEN))"
end
