# ----------------------------
# Warping safety parameters
# ----------------------------

"""
    CaseParams

Safety parameters for a specific warping case (NNZZ or NZZZ).

# Fields
- `threshold_distance::Float64`: Maximum distance of centroid from isosurface for warping
- `max_node_displacement::Float64`: Maximum allowed node displacement during warping
- `min_volume_ratio::Float64`: Minimum volume ratio relative to original element
- `min_dihedral_angle::Float64`: Minimum allowed dihedral angle in degrees
- `max_dihedral_angle::Float64`: Maximum allowed dihedral angle in degrees
"""
struct CaseParams
    threshold_distance::Float64
    max_node_displacement::Float64
    min_volume_ratio::Float64
    min_dihedral_angle::Float64
    max_dihedral_angle::Float64
end

"""
    WarpingSafetyParams

Container for safety parameters for different warping cases.

# Fields
- `nnzz::CaseParams`: Parameters for NNZZ case (two nodes outside, two on surface)
- `nzzz::CaseParams`: Parameters for NZZZ case (one node outside, three on surface)
"""
struct WarpingSafetyParams
    nnzz::CaseParams
    nzzz::CaseParams
end

"""
    create_warping_params(scheme::String, grid_step::Float64) -> WarpingSafetyParams

Create scheme-specific safety parameters for warping operations.

NNZZ case (2 nodes warped) uses stricter parameters.
NZZZ case (1 node warped) uses relaxed parameters.

# Arguments
- `scheme::String`: Discretization scheme ("A15" or "Schlafli")
- `grid_step::Float64`: Mesh grid step size

# Returns
- `WarpingSafetyParams`: Safety parameters for warping operations
"""
function create_warping_params(scheme::String, grid_step::Float64)
    if scheme == "A15"
        base_threshold = 0.15 * grid_step
        # NNZZ: Strict parameters (warping 2 nodes simultaneously)
        nnzz = CaseParams(base_threshold, 0.2 * grid_step, 0.05, 10.0, 140.0)
        # NZZZ: Relaxed parameters (warping only 1 node)
        nzzz = CaseParams(2.0 * base_threshold, 0.4 * grid_step, 0.025, 10.0, 140.0)
    elseif scheme == "Schlafli"
        base_threshold = 0.3 * grid_step
        nnzz = CaseParams(base_threshold, 0.4 * grid_step, 0.05, 10.0, 140.0)
        nzzz = CaseParams(2.0 * base_threshold, 0.8 * grid_step, 0.025, 10.0, 140.0)
    else
        error("Unknown scheme: $scheme")
    end
    return WarpingSafetyParams(nnzz, nzzz)
end

"""
    compute_dihedral_angle_range(mesh::BlockMesh, tet::Vector{Int}) -> (Float64, Float64)

Compute the minimum and maximum interior dihedral angles for a tetrahedron.
Dihedral angle is measured between two adjacent faces of the tetrahedron.

# Returns
Tuple of (min_angle, max_angle) in degrees
"""
function compute_dihedral_angle_range(mesh::BlockMesh, tet::Vector{Int})
    # Get tetrahedron vertices
    vertices = [mesh.X[tet[i]] for i = 1:4]

    # Compute outward-facing face normals for all 4 faces
    face_normals = [
        normalize(cross(vertices[2] - vertices[1], vertices[3] - vertices[1])),  # Face 1-2-3
        normalize(cross(vertices[2] - vertices[1], vertices[4] - vertices[1])),  # Face 1-2-4
        normalize(cross(vertices[3] - vertices[1], vertices[4] - vertices[1])),  # Face 1-3-4
        normalize(cross(vertices[3] - vertices[2], vertices[4] - vertices[2])),  # Face 2-3-4
    ]

    # Ensure normals point outward from tetrahedron
    opposite_vertices = [vertices[4], vertices[3], vertices[2], vertices[1]]
    for i = 1:4
        # Calculate face center
        face_indices = setdiff(1:4, [i])
        face_center = sum(vertices[j] for j in face_indices) / 3.0

        # Vector from face center to opposite vertex
        to_opposite = opposite_vertices[i] - face_center

        # Flip normal if it points inward
        if dot(face_normals[i], to_opposite) > 0
            face_normals[i] = -face_normals[i]
        end
    end

    # Compute all 6 interior dihedral angles
    # Interior angle = π - angle_between_outward_normals = acos(-dot(n1, n2))
    min_angle = 180.0  # Initialize to maximum possible
    max_angle = 0.0    # Initialize to minimum possible

    for i = 1:3
        for j = (i+1):4
            # Compute interior dihedral angle
            interior_angle = acos(-clamp(dot(face_normals[i], face_normals[j]), -1.0, 1.0))
            # Convert to degrees
            angle_deg = interior_angle * 180.0 / π

            # Update min and max
            min_angle = min(min_angle, angle_deg)
            max_angle = max(max_angle, angle_deg)
        end
    end

    return (min_angle, max_angle)
end
