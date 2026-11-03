# ----------------------------------------------------------------------------
# Signed distance of the design domain (domain clip of rho2sdf.jl)
# ----------------------------------------------------------------------------
# A port of rho2sdf.jl's src/SdfSmoothing/DomainDistance.jl. rho2sdf.jl sets clip_shift so that
# its clipped body keeps the volume; the mesher must evaluate the SAME d_D for the meshed body to
# keep it, so the two files differ only in the input (the matrices of the JLD2 instead of
# rho2sdf.jl's `Mesh`). Keep them in step.
#
# Signed distance d_D of the design domain D, positive inside as the SDF of rho2sdf.
#
# D is the region covered by the hexahedral mesh of the SIMP result. The RBF field of the body is
# clipped by it, φ = min(φ̃ - shift, d_D) (`RBFs_smoothing`), so that the body stays inside D.
#
# The boundary of D consists of the element faces that belong to one element only. Each face is
# split into 4 triangles around its centre (`face_triangles`):
# - |d_D| is the distance to the nearest boundary triangle;
# - d_D > 0 if the point lies in an element, each element being split into 24 tetrahedra between
#   its centre and the triangles of its faces.
# The distance and the sign thus describe the same polyhedron. For planar faces (a box) it is the
# mesh itself and d_D is exact; warped faces are approximated by their 4 triangles.

const Point = SVector{3,Float64}

# Local nodes of the 6 faces of a HEX8 element in VTK corner order (rho2sdf.jl `Mesh.ISN`)
const HEX8_FACES = [[1, 4, 3, 2], [1, 2, 6, 5], [2, 3, 7, 6], [3, 4, 8, 7], [4, 1, 5, 8], [5, 6, 7, 8]]

"""
    face_triangles(corners, face)

The 4 triangles (corner, next corner, face centre) of the quadrilateral face `face` (4 local
node numbers, `HEX8_FACES`) of an element with corners `corners`.
"""
function face_triangles(corners::NTuple{8,Point}, face::Vector{Int})
    a, b, c, d = corners[face[1]], corners[face[2]], corners[face[3]], corners[face[4]]
    centre = (a + b + c + d) / 4
    return ((a, b, centre), (b, c, centre), (c, d, centre), (d, a, centre))
end

"""
    closest_point_on_triangle(x, a, b, c)

Point of the triangle (a, b, c) nearest to x. The regions of the vertices, the edges and the
interior are tested in turn (Ericson, Real-Time Collision Detection, sec. 5.1.5).
"""
function closest_point_on_triangle(x::Point, a::Point, b::Point, c::Point)
    ab = b - a
    ac = c - a

    # Vertex a:
    d1 = dot(ab, x - a)
    d2 = dot(ac, x - a)
    if d1 <= 0 && d2 <= 0
        return a
    end

    # Vertex b:
    d3 = dot(ab, x - b)
    d4 = dot(ac, x - b)
    if d3 >= 0 && d4 <= d3
        return b
    end

    # Edge ab:
    vc = d1 * d4 - d3 * d2
    if vc <= 0 && d1 >= 0 && d3 <= 0
        return a + d1 / (d1 - d3) * ab
    end

    # Vertex c:
    d5 = dot(ab, x - c)
    d6 = dot(ac, x - c)
    if d6 >= 0 && d5 <= d6
        return c
    end

    # Edge ac:
    vb = d5 * d2 - d1 * d6
    if vb <= 0 && d2 >= 0 && d6 <= 0
        return a + d2 / (d2 - d6) * ac
    end

    # Edge bc:
    va = d3 * d6 - d5 * d4
    if va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0
        return b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b)
    end

    # Interior:
    return a + (vb * ab + vc * ac) / (va + vb + vc)
end

"""Six times the signed volume of the tetrahedron (a, b, c, d)."""
orientation(a::Point, b::Point, c::Point, d::Point) = dot(b - a, cross(c - a, d - a))

"""
    in_tetrahedron(x, a, b, c, d)

Whether x lies in the tetrahedron (a, b, c, d): its barycentric coordinates are all ≥ 0. A small
tolerance keeps points on a face shared by two tetrahedra inside despite round-off.
"""
function in_tetrahedron(x::Point, a::Point, b::Point, c::Point, d::Point)
    tolerance = 1e-12
    volume = orientation(a, b, c, d)
    if orientation(x, b, c, d) / volume < -tolerance
        return false
    end
    if orientation(a, x, c, d) / volume < -tolerance
        return false
    end
    if orientation(a, b, x, d) / volume < -tolerance
        return false
    end
    return orientation(a, b, c, x) / volume >= -tolerance
end

"""
    in_element(x, corners, faces)

Whether x lies in the hexahedral element with corners `corners`, split into 24 tetrahedra
between its centre and the triangles of its faces (`face_triangles`).
"""
function in_element(x::Point, corners::NTuple{8,Point}, faces::Vector{Vector{Int}})
    centre = sum(corners) / 8
    for face in faces
        for (a, b, face_centre) in face_triangles(corners, face)
            if in_tetrahedron(x, centre, a, b, face_centre)
                return true
            end
        end
    end
    return false
end

"""
    DomainDistance(nodes, elements, cap)

Signed distance d_D of the design domain covered by a hexahedral mesh, positive inside (see the
top of the file). Usage: `d_D = DomainDistance(nodes, elements, cap)`, then `d_D(x)`.

1. Boundary faces: the element faces that belong to one element only, 4 triangles each.
2. Bins: cubes of edge cap/2 over the box of the mesh enlarged by `cap`. Each bin lists the
   triangles that can lie within `cap` of it, sorted by their distance from the bin centre, and
   the elements whose bounding box overlaps it.
3. `d_D(x)`: distance to the nearest triangle of the bin of x, at most `cap`; positive if x lies
   in one of the elements of the bin (`in_element`).

|d_D| is exact up to `cap`. Farther points get ±`cap`: only their sign is used.

# Arguments
- `nodes::AbstractMatrix{<:Real}`: Node coordinates, 3 × nodes (rho2sdf.jl key `domain_nodes`)
- `elements::AbstractMatrix{<:Integer}`: HEX8 elements in VTK corner order, 8 × elements
  (rho2sdf.jl key `domain_elements`)
- `cap::Real`: Largest |d_D| returned
"""
struct DomainDistance
    corners::Vector{NTuple{8,Point}}           # corners of each element
    element_boxes::Vector{NTuple{2,Point}}     # bounding box (min, max) of each element
    faces::Vector{Vector{Int}}                 # local nodes of the element faces (HEX8_FACES)
    triangles::Vector{NTuple{3,Point}}         # boundary faces, 4 triangles each
    cap::Float64                               # largest |d_D| returned
    bin_origin::Point                          # lower corner of the first bin
    bin_size::Float64                          # edge of a bin
    bin_dims::NTuple{3,Int}                    # number of bins along x, y, z
    bin_triangles::Vector{Vector{Int}}         # triangles that can lie within cap of each bin
    bin_centre_distances::Vector{Vector{Float64}}  # their distances from the bin centre, ascending
    bin_elements::Vector{Vector{Int}}          # elements whose bounding box overlaps each bin
end

function DomainDistance(nodes::AbstractMatrix{<:Real}, elements::AbstractMatrix{<:Integer}, cap::Real)
    size(nodes, 1) == 3 || error("DomainDistance: nodes must be 3 × nodes, got $(size(nodes))")
    size(elements, 1) == 8 || error("DomainDistance: elements must be 8 × elements, got $(size(elements))")
    all(1 .<= elements .<= size(nodes, 2)) ||
        error("DomainDistance: element node numbers outside 1:$(size(nodes, 2))")
    cap > 0 || error("DomainDistance: cap must be positive, got $cap")
    cap = Float64(cap)
    n_elements = size(elements, 2)

    points = [Point(nodes[1, n], nodes[2, n], nodes[3, n]) for n = 1:size(nodes, 2)]
    corners = [ntuple(a -> points[elements[a, e]], 8) for e = 1:n_elements]
    element_boxes = [(min.(c...), max.(c...)) for c in corners]

    # 1. Boundary faces, found by counting the faces by their sorted node numbers:
    face_key(e, face) = Tuple(sort(SVector{4,Int}(elements[face, e])))
    face_count = Dict{NTuple{4,Int},Int}()
    for e = 1:n_elements
        for face in HEX8_FACES
            key = face_key(e, face)
            face_count[key] = get(face_count, key, 0) + 1
        end
    end
    triangles = NTuple{3,Point}[]
    for e = 1:n_elements
        for face in HEX8_FACES
            if face_count[face_key(e, face)] == 1
                append!(triangles, face_triangles(corners[e], face))
            end
        end
    end

    # 2. Bins over the box of the mesh enlarged by cap:
    bin_size = cap / 2
    bin_origin = reduce((p, q) -> min.(p, q), points) .- cap
    bin_end = reduce((p, q) -> max.(p, q), points) .+ cap
    bin_dims = Tuple(ceil.(Int, (bin_end .- bin_origin) ./ bin_size))
    n_bins = prod(bin_dims)
    bin_triangles = [Int[] for _ = 1:n_bins]
    bin_centre_distances = [Float64[] for _ = 1:n_bins]
    bin_elements = [Int[] for _ = 1:n_bins]
    first_bin(x) = max.(1, floor.(Int, (x .- bin_origin) ./ bin_size) .+ 1)
    last_bin(x) = min.(bin_dims, floor.(Int, (x .- bin_origin) ./ bin_size) .+ 1)
    linear_bin(i, j, k) = i + bin_dims[1] * (j - 1 + bin_dims[2] * (k - 1))

    # A triangle can lie within cap of a bin if it lies within cap + half the bin diagonal of its centre:
    reach = cap + sqrt(3) / 2 * bin_size
    for (t, (a, b, c)) in enumerate(triangles)
        low, high = first_bin(min.(a, b, c) .- cap), last_bin(max.(a, b, c) .+ cap)
        for k = low[3]:high[3]
            for j = low[2]:high[2]
                for i = low[1]:high[1]
                    bin_centre = bin_origin + bin_size * Point(i - 0.5, j - 0.5, k - 0.5)
                    centre_distance = norm(bin_centre - closest_point_on_triangle(bin_centre, a, b, c))
                    if centre_distance <= reach
                        push!(bin_triangles[linear_bin(i, j, k)], t)
                        push!(bin_centre_distances[linear_bin(i, j, k)], centre_distance)
                    end
                end
            end
        end
    end
    for bin = 1:n_bins
        order = sortperm(bin_centre_distances[bin])
        bin_triangles[bin] = bin_triangles[bin][order]
        bin_centre_distances[bin] = bin_centre_distances[bin][order]
    end

    for (e, (box_min, box_max)) in enumerate(element_boxes)
        low, high = first_bin(box_min), last_bin(box_max)
        for k = low[3]:high[3]
            for j = low[2]:high[2]
                for i = low[1]:high[1]
                    push!(bin_elements[linear_bin(i, j, k)], e)
                end
            end
        end
    end

    return DomainDistance(corners, element_boxes, HEX8_FACES, triangles, cap,
        bin_origin, bin_size, bin_dims, bin_triangles, bin_centre_distances, bin_elements)
end

"""
    (d_D::DomainDistance)(x) -> Float64

Signed distance of the design domain at x, positive inside, |d_D| ≤ cap.
"""
function (d_D::DomainDistance)(x::Point)
    ijk = floor.(Int, (x .- d_D.bin_origin) ./ d_D.bin_size) .+ 1
    if any(ijk .< 1) || any(ijk .> d_D.bin_dims)
        return -d_D.cap     # outside the bins: farther than cap from the mesh
    end
    bin = ijk[1] + d_D.bin_dims[1] * (ijk[2] - 1 + d_D.bin_dims[2] * (ijk[3] - 1))

    # 1. Distance to the nearest boundary triangle. A triangle is at least its distance from the bin
    #    centre minus |x - centre| away from x, so the sorted search stops at the first one beyond:
    bin_centre = d_D.bin_origin + d_D.bin_size * (ijk .- 0.5)
    centre_offset = norm(x - bin_centre)
    distance = d_D.cap
    for (t, centre_distance) in zip(d_D.bin_triangles[bin], d_D.bin_centre_distances[bin])
        if centre_distance - centre_offset >= distance
            break
        end
        a, b, c = d_D.triangles[t]
        distance = min(distance, norm(x - closest_point_on_triangle(x, a, b, c)))
    end

    # 2. Sign: inside if x lies in one of the elements of the bin:
    for e in d_D.bin_elements[bin]
        box_min, box_max = d_D.element_boxes[e]
        if all(box_min .<= x .<= box_max) && in_element(x, d_D.corners[e], d_D.faces)
            return distance
        end
    end
    return -distance
end
