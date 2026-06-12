# ----------------------------
# Surface-tet acceptance parameters
# ----------------------------

"""
    DihedralBounds

Acceptable interior dihedral-angle range (in degrees) for a retained surface
("quadruple-zero") tetrahedron. A candidate whose minimum interior angle is below
`min_dihedral_angle`, or whose maximum interior angle is above `max_dihedral_angle`,
is too poorly shaped and gets discarded (Labelle §3.4 surface-tet quality test).

# Fields
- `min_dihedral_angle::Float64`: Minimum allowed dihedral angle in degrees
- `max_dihedral_angle::Float64`: Maximum allowed dihedral angle in degrees
"""
struct DihedralBounds
    min_dihedral_angle::Float64
    max_dihedral_angle::Float64
end

"""
    create_warping_params(scheme::String) -> DihedralBounds

Return the scheme-specific dihedral-angle bounds used when deciding whether to keep a
surface ("quadruple-zero") tetrahedron during slicing. Only "A15" is supported.

# Arguments
- `scheme::String`: Discretization scheme (only "A15" is supported)

# Returns
- `DihedralBounds`: the acceptable interior dihedral-angle range
"""
function create_warping_params(scheme::String)
    if scheme != "A15"
        error("Unknown scheme: $scheme. Only 'A15' is supported.")
    end
    return DihedralBounds(10.0, 140.0)
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
