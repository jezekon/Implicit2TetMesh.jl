# Functions for tetrahedral mesh slicing along an isosurface

"""
    slice_ambiguous_tetrahedra!(mesh::BlockMesh, scheme::String; cut_points = :linear)

Slice tetrahedra that stick out across the isosurface (SDF zero level set), trimming each one
back to the interior region (SDF <= 0). This is a direct port of quartet's `trim_spikes`
(make_tet_mesh.cpp:221-346). Every spike tetrahedron is replaced by smaller tetrahedra whose
vertices lie either on the original lattice or on the edge cut points produced by `cut_edge!`.

`cut_points` selects how the cut point on a crossing edge is located (`:linear` = quartet's
endpoint-value estimate, the default; `:bisection` = the true zero of `eval_sdf` along the
edge, Labelle & Shewchuk §3.1 -- required for input fields that are not distance-like). It
MUST match the mode `warp!` was run with, otherwise the warped vertices and the cut vertices
disagree about where the surface is.

The "quadruple-zero" tetrahedra (all four vertices on the surface) are handled separately in a
SECOND pass by `resolve_surface_candidates!`, following Labelle's surface-fidelity heuristic
(§3.4). They have to be deferred because their fate depends on whether their faces adjoin the rest
of the solid mesh, which is only known once the whole solid mesh has been built.

# Arguments
- `mesh::BlockMesh`: The mesh to process
- `scheme::String`: Discretization scheme (only "A15" is supported)
"""
function slice_ambiguous_tetrahedra!(
    mesh::BlockMesh,
    scheme::String;
    cut_points::Symbol = :linear,
)
    @info "Slicing tetrahedra using trim_spikes logic (cut_points = $cut_points)..."
    validate_cut_points(cut_points)

    cut_map = Dict{Tuple{Int,Int},Int}()
    solid_tets = Vector{Vector{Int64}}()          # definitely-interior output tets (no quadruple-zero)
    surface_candidates = Vector{Vector{Int64}}()  # all-on-surface tets, resolved in the second pass
    sizehint!(solid_tets, length(mesh.IEN))

    current_IEN = mesh.IEN
    mesh.IEN = Vector{Vector{Int64}}()

    # First pass: trim every tetrahedron. Crossing and fully-interior tets are emitted straight
    # into the solid mesh; all-on-surface ("quadruple-zero") tets are collected as candidates.
    for tet in current_IEN
        apply_stencil_trim_spikes!(mesh, tet, cut_map, solid_tets, surface_candidates, cut_points)
    end

    # Second pass: decide which deferred surface-only tets to keep (Labelle §3.4 heuristic).
    retained = resolve_surface_candidates!(mesh, solid_tets, surface_candidates, scheme)

    mesh.IEN = solid_tets
    append!(mesh.IEN, retained)
    println(
        "  After slicing: $(length(mesh.IEN)) tetrahedra " *
        "($(length(retained))/$(length(surface_candidates)) surface tets retained)",
    )
end

"""
    cut_edge!(i, j, mesh, node_sdf, cut_map, cut_points = :linear) -> Int

Find or create the vertex where the isosurface crosses the edge `(i, j)`, placed at the cut
point

    X_cut = (1 - alpha) * X_i + alpha * X_j,    alpha in (0, 1).

With `cut_points = :linear` (default) `alpha = phi_i / (phi_i - phi_j)`, exactly as in
quartet's `cut_edge` (make_tet_mesh.cpp:196-217) -- exact for distance-like fields. With
`:bisection` the true zero of `eval_sdf` along the edge is used (`bisect_cut_alpha`,
Labelle & Shewchuk §3.1) -- required for non-distance inputs, where the linear estimate puts
the cut vertex visibly inside the solid. Either way the vertex stays ON the edge (no
projection off it), which keeps the trimmed tetrahedra well shaped, and the mode must be
consistent with the one used by the edge-based `warp!`. Each edge is cut only once: the new
vertex index is cached in `cut_map` (keyed by the canonical edge), so all tetrahedra sharing
the edge reuse the same vertex and the boundary stays watertight.

If an endpoint already lies on the surface (`|phi| < tol`) that endpoint index is returned
instead of creating a new vertex.
"""
function cut_edge!(
    i::Int,
    j::Int,
    mesh::BlockMesh,
    node_sdf::Vector{Float64},
    cut_map::Dict{Tuple{Int,Int},Int},
    cut_points::Symbol = :linear,
)::Int
    # Get SDF values for both endpoints
    sdf_i = node_sdf[i]
    sdf_j = node_sdf[j]

    tol = mesh.grid_tol
    is_i_zero = abs(sdf_i) < tol
    is_j_zero = abs(sdf_j) < tol

    # Handle cases where an endpoint is already on the surface
    if is_i_zero && is_j_zero
        return min(i, j)
    elseif is_i_zero
        return i
    elseif is_j_zero
        return j
    end

    # The two endpoints must lie on opposite sides of the surface
    if sign(sdf_i) == sign(sdf_j)
        error(
            "cut_edge! called with nodes on same side: i=$i (sdf=$sdf_i), j=$j (sdf=$sdf_j)",
        )
    end

    # Canonical edge representation so the edge is cut only once
    edge = (min(i, j), max(i, j))
    if haskey(cut_map, edge)
        return cut_map[edge]
    end

    # Cut point along the edge: alpha is the fraction from i to the zero crossing, and the
    # point is the straight-line interpolation between the endpoints. :linear estimates alpha
    # from the endpoint values (quartet cut_edge); :bisection locates the actual zero of
    # eval_sdf on the segment.
    alpha =
        cut_points === :bisection ?
        bisect_cut_alpha(mesh, mesh.X[i], mesh.X[j], sdf_i, sdf_j) :
        sdf_i / (sdf_i - sdf_j)
    cut_point = (1.0 - alpha) * mesh.X[i] + alpha * mesh.X[j]

    push!(mesh.X, cut_point)
    push!(mesh.node_sdf, 0.0)   # the cut point lies exactly on the surface
    new_index = length(mesh.X)
    cut_map[edge] = new_index
    return new_index
end

"""
    is_positively_oriented(mesh::BlockMesh, tet) -> Bool

The single, EXACT orientation test for a tetrahedron. Returns true iff the tetrahedron
`(v1, v2, v3, v4)` has strictly positive signed volume, i.e. the Jacobian determinant
`det[v2-v1, v3-v1, v4-v1] > 0` -- the same sign convention the code used before, but decided
EXACTLY (no floating-point tolerance) via Shewchuk's adaptive predicate `ExactPredicates.orient`.

NOTE on the sign: `ExactPredicates.orient(a, b, c, d)` returns the OPPOSITE sign of
`dot(b-a, cross(c-a, d-a))` (verified empirically), so a positive signed volume corresponds to
`orient(...) < 0`. A flat (zero-volume / coplanar) tetrahedron returns false.

This is the ONE place the orientation sign is decided; every orientation/inversion check in this
file routes through it so the decision is never duplicated.
"""
function is_positively_oriented(mesh::BlockMesh, tet)::Bool
    return orient(mesh.X[tet[1]], mesh.X[tet[2]], mesh.X[tet[3]], mesh.X[tet[4]]) < 0
end

"""
    four_distinct(a, b, c, d) -> Bool

True iff the four node indices are pairwise distinct, i.e. the tetrahedron is non-degenerate.
Equivalent to `length(Set((a, b, c, d))) == 4` but without allocating a `Set` -- this runs once
per emitted/checked tetrahedron, so it is kept allocation-free.
"""
function four_distinct(a::Int, b::Int, c::Int, d::Int)::Bool
    return a != b && a != c && a != d && b != c && b != d && c != d
end

"""
    vert_greater(vert_data, i, j) -> Bool

Strict ordering used to sort a crossing tetrahedron's vertices, matching quartet `trim_spikes`:
descending by SDF, with ties broken so the LARGER original node index is treated as the
more-outside vertex (the same direction quartet uses when it swaps on `index <`). Because this is
a global total order on (SDF, index), face-adjacent tetrahedra sort a shared quad identically and
therefore split it along the same diagonal -- this is what keeps the boundary watertight.

`vert_data` is the tuple of `(sdf, node-index)` pairs; passing it as an argument keeps this a plain
function so the sort below needs no per-tetrahedron closure allocation.
"""
function vert_greater(vert_data, i::Int, j::Int)::Bool
    return (vert_data[i][1] > vert_data[j][1]) ||
           (vert_data[i][1] == vert_data[j][1] && vert_data[i][2] > vert_data[j][2])
end

"""
    check_tetrahedron_orientation(mesh::BlockMesh, tet::Vector{Int})

Check if a tetrahedron has positive orientation (positive Jacobian determinant).
Returns true for correctly oriented tetrahedra, false for inverted (or flat) ones. Validates the
indices first, then delegates the sign decision to the exact `is_positively_oriented`.
"""
function check_tetrahedron_orientation(mesh::BlockMesh, tet::Vector{Int})
    # Validate indices
    if length(tet) != 4 || any(i -> i <= 0 || i > length(mesh.X), tet)
        @warn "Invalid tetrahedron indices for orientation check: $tet"
        return false
    end

    return is_positively_oriented(mesh, tet)
end

"""
    fix_tetrahedron_orientation!(mesh::BlockMesh, tet::Vector{Int})

Fix the orientation of an inverted tetrahedron by swapping vertices.
Modifies the tetrahedron in-place. Returns true if fixed, false otherwise.
"""
function fix_tetrahedron_orientation!(mesh::BlockMesh, tet::Vector{Int})
    if length(tet) != 4
        @warn "Attempting to fix orientation of non-tetrahedron: $tet"
        return false
    end

    if !check_tetrahedron_orientation(mesh, tet)
        # Swap last two vertices to flip orientation
        tet[3], tet[4] = tet[4], tet[3]
        return true  # Orientation was fixed
    end

    return false # No fix needed
end

"""
    apply_stencil_trim_spikes!(mesh, tet, cut_map, solid_tets, surface_candidates,
                               cut_points = :linear)

Apply the trimming algorithm to a single tetrahedron that crosses the isosurface.
The trimmed interior tetrahedra (SDF <= 0) are pushed straight into `solid_tets`; there is no
per-tetrahedron return array (this function runs once for every tetrahedron in the mesh, so
avoiding that allocation matters). The SDF sign pattern is classified case-by-case, ported 1:1
from quartet `trim_spikes` (make_tet_mesh.cpp:221-346). `cut_points` is forwarded to
`cut_edge!` (see there).

A "quadruple-zero" tetrahedron (all four vertices on the surface) is NOT decided here: it is
appended to `surface_candidates` and resolved later by `resolve_surface_candidates!` (Labelle
§3.4), because its fate depends on the surrounding solid mesh.
"""
function apply_stencil_trim_spikes!(
    mesh::BlockMesh,
    tet::Vector{Int64},
    cut_map::Dict{Tuple{Int,Int},Int},
    solid_tets::Vector{Vector{Int64}},
    surface_candidates::Vector{Vector{Int64}},
    cut_points::Symbol = :linear,
)
    # SDF value at each of the four nodes. Read as scalars (no per-tet temporary array): the fast
    # paths below run for EVERY tetrahedron, and the vast majority are fully interior, so a
    # temporary here would dominate the slicing cost.
    tol = mesh.grid_tol
    s1 = mesh.node_sdf[tet[1]]
    s2 = mesh.node_sdf[tet[2]]
    s3 = mesh.node_sdf[tet[3]]
    s4 = mesh.node_sdf[tet[4]]

    # --- Special cases handling ---
    # Case 1: All nodes outside (SDF > tol) - discard tetrahedron (all-outside case)
    if s1 > tol && s2 > tol && s3 > tol && s4 > tol
        return
    end

    # Case 2: All nodes inside (SDF < -tol) - keep original (all-inside case)
    if s1 < -tol && s2 < -tol && s3 < -tol && s4 < -tol
        push!(solid_tets, tet)
        return
    end

    # Case 3: All nodes on the surface (|SDF| <= tol) - "quadruple-zero" tetrahedron. Its fate is
    # ambiguous (Labelle §3.4) and depends on the surrounding solid mesh, so defer it: collect it
    # as a candidate and emit nothing now. resolve_surface_candidates! decides it afterwards.
    if abs(s1) <= tol && abs(s2) <= tol && abs(s3) <= tol && abs(s4) <= tol
        push!(surface_candidates, tet)
        return
    end

    # --- General case: Tetrahedron crossing the isosurface ---
    # Sort vertices by SDF value (largest first, i.e. most-outside vertex first). vert_data is a
    # tuple of (sdf, node-index) pairs, so the sort below indexes it without allocating.
    vert_data = ((s1, tet[1]), (s2, tet[2]), (s3, tet[3]), (s4, tet[4]))
    p = MVector(1, 2, 3, 4)   # mutable permutation, stack-allocated (no per-tet heap allocation)
    flipped = false

    # Sort the four vertices into descending (SDF, index) order with vert_greater, tracking
    # orientation flips. (See vert_greater for the tie-break that keeps the boundary watertight.)
    if vert_greater(vert_data, p[2], p[1])
        p[1], p[2] = p[2], p[1]
        flipped = !flipped
    end
    if vert_greater(vert_data, p[4], p[3])
        p[3], p[4] = p[4], p[3]
        flipped = !flipped
    end
    if vert_greater(vert_data, p[3], p[1])
        p[1], p[3] = p[3], p[1]
        flipped = !flipped
    end
    if vert_greater(vert_data, p[4], p[2])
        p[2], p[4] = p[4], p[2]
        flipped = !flipped
    end
    if vert_greater(vert_data, p[3], p[2])
        p[2], p[3] = p[3], p[2]
        flipped = !flipped
    end

    # Sorted vertices (s has highest SDF = most outside, p has lowest = most inside).
    # Index the tuple directly instead of building temporary arrays.
    s_idx = vert_data[p[1]][2]
    r_idx = vert_data[p[2]][2]
    q_idx = vert_data[p[3]][2]
    p_idx = vert_data[p[4]][2]
    sdf_s = vert_data[p[1]][1]
    sdf_r = vert_data[p[2]][1]
    sdf_q = vert_data[p[3]][1]
    sdf_p = vert_data[p[4]][1]

    # Check if SDF values are properly sorted (descending)
    if !(sdf_s >= sdf_r >= sdf_q >= sdf_p)
        @warn "SDF values not sorted: s=$sdf_s, r=$sdf_r, q=$sdf_q, p=$sdf_p"
    end

    # (All-on-surface "quadruple-zero" tets were already deferred to surface_candidates above,
    #  so by here at least one vertex is strictly inside or strictly outside.)

    # Classify each vertex by geometric role under the phi < 0 = inside convention:
    #   outside : sdf >  tol   (positive)
    #   inside  : sdf < -tol   (negative)
    #   surface : |sdf| <= tol (zero)
    # The case labels below (NNNP, NPPP, ...) name the geometric pattern and line up 1:1 with
    # quartet's notation (+++-, +---, ...); see the per-branch comments.
    is_s_outside = sdf_s > tol
    is_r_outside = sdf_r > tol
    is_q_outside = sdf_q > tol

    is_p_inside = sdf_p < -tol
    is_q_inside = sdf_q < -tol
    is_r_inside = sdf_r < -tol

    is_s_surface = abs(sdf_s) <= tol
    is_r_surface = abs(sdf_r) <= tol
    is_q_surface = abs(sdf_q) <= tol
    is_p_surface = abs(sdf_p) <= tol

    # Helper: emit a trimmed tetrahedron straight into the solid mesh, with an orientation check.
    function add_tet!(t::Vector{Int})
        # Skip degenerate cases (duplicate vertices)
        if !four_distinct(t[1], t[2], t[3], t[4])
            return
        end

        # Fix orientation and emit
        fix_tetrahedron_orientation!(mesh, t)
        if !check_tetrahedron_orientation(mesh, t)
            @warn "Tetrahedron $t still has incorrect orientation after fix. Skipping."
            return
        end
        push!(solid_tets, t)
    end

    # --- Case analysis based on SDF sign patterns --- (S >= R >= Q >= P)

    # Entirely outside or on the surface (quartet's vphi[s]==0 branch: +++0 / ++00 / +000).
    # The most-inside vertex is already on the surface, so no interior region remains -> discard.
    if is_p_surface
        return

        # Surface/interior tetrahedron with no outside vertex: ZZZP, ZZPP, ZPPP. quartet leaves
        # these untouched (all phi <= 0), so we keep them as-is.
    elseif is_s_surface && !is_r_outside && !is_q_outside && is_p_inside
        push!(solid_tets, tet)
        return

        # Case NNNP (quartet +++-): three nodes outside, one inside
    elseif is_s_outside && is_r_outside && is_q_outside && is_p_inside
        # Cut three edges from outside nodes to inside node
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)
        qp = cut_edge!(q_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)

        # Create one tetrahedron from three cut points and interior node
        if flipped
            add_tet!([rp, sp, qp, p_idx])
        else
            add_tet!([sp, rp, qp, p_idx])
        end

        # Case NPPP (quartet +---): one node outside, three inside
    elseif is_s_outside && is_r_inside
        if !is_r_outside && !is_q_outside # Confirm r, q, p are interior nodes
            # Cut edges from outside node to each interior node
            sr = cut_edge!(s_idx, r_idx, mesh, mesh.node_sdf, cut_map, cut_points)
            sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map, cut_points)
            sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)

            # Tetrahedralize the resulting triangular prism. The quad faces are split to the
            # deepest vertex; the consistent sort above makes the split match face-adjacent tets.
            if flipped
                add_tet!([q_idx, r_idx, p_idx, sr])
                add_tet!([p_idx, q_idx, sr, sq])
                add_tet!([sr, p_idx, sq, sp])
            else
                add_tet!([r_idx, q_idx, p_idx, sr])
                add_tet!([q_idx, p_idx, sr, sq])
                add_tet!([p_idx, sr, sq, sp])
            end
        else
            @warn "Logic error in NPPP branch: r=$sdf_r, q=$sdf_q, p=$sdf_p"
            return
        end

        # Case NNPP (quartet ++--): two nodes outside, two inside
    elseif is_r_outside && is_q_inside
        # Cut all four edges crossing the isosurface
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map, cut_points)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)
        rq = cut_edge!(r_idx, q_idx, mesh, mesh.node_sdf, cut_map, cut_points)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)

        # Create three tetrahedra for interior region (quad split consistent with the sort)
        if flipped
            add_tet!([p_idx, rq, q_idx, sq])
            add_tet!([p_idx, sp, sq, rp])
            add_tet!([p_idx, rp, sq, rq])
        else
            add_tet!([p_idx, q_idx, rq, sq])
            add_tet!([p_idx, sq, sp, rp])
            add_tet!([p_idx, sq, rp, rq])
        end

        # Case NZPP (quartet +0--): one outside, one on surface, two inside
    elseif is_s_outside && is_r_surface && is_q_inside
        # Cut edges from outside node to interior nodes
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map, cut_points)

        # Create two tetrahedra
        if flipped
            add_tet!([q_idx, r_idx, p_idx, sq])
            add_tet!([p_idx, r_idx, sq, sp])
        else
            add_tet!([r_idx, q_idx, p_idx, sq])
            add_tet!([r_idx, p_idx, sq, sp])
        end

        # Case NNZP (quartet ++0-): two outside, one on surface, one inside
    elseif is_r_outside && is_q_surface && is_p_inside
        # Cut edges from outside nodes to inside node
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)

        # Create one tetrahedron
        if flipped
            add_tet!([p_idx, q_idx, rp, sp])
        else
            add_tet!([q_idx, p_idx, rp, sp])
        end

        # Case NZZP (quartet +00-): one outside, two on surface, one inside
    elseif is_s_outside && is_r_surface && is_q_surface && is_p_inside
        # Cut edge from outside node to inside node
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map, cut_points)

        # Create one tetrahedron
        if flipped
            add_tet!([q_idx, r_idx, p_idx, sp])
        else
            add_tet!([r_idx, q_idx, p_idx, sp])
        end

    else
        @warn "Unexpected SDF pattern in apply_stencil_trim_spikes!: s=$sdf_s, r=$sdf_r, q=$sdf_q, p=$sdf_p. Indices: s=$s_idx, r=$r_idx, q=$q_idx, p=$p_idx. Tet: $tet. Flipped: $flipped"
        return # Discard in case of unexpected pattern
    end

    return
end

"""
    face_key(a::Int, b::Int, c::Int) -> NTuple{3,Int}

Return the three node indices sorted ascending, as a tuple, WITHOUT allocating a temporary array
(a small 3-element sorting network). This is the canonical orientation-independent key for a
triangular face. It mirrors the identical helper in `src/Modification/RemoveIsolatedComponents.jl`;
the function is duplicated here because `GenerateMesh` is included before `Modification` and so
cannot depend on it.
"""
function face_key(a::Int, b::Int, c::Int)::NTuple{3,Int}
    a > b && ((a, b) = (b, a))
    b > c && ((b, c) = (c, b))
    a > b && ((a, b) = (b, a))
    return (a, b, c)
end

"""
    tetrahedron_faces(tet) -> NTuple{4, NTuple{3,Int}}

Return the four triangular faces of a tetrahedron as sorted node-index triples. Sorting (via
`face_key`) makes each key orientation-independent, so a face shared by two tetrahedra maps to the
same key (the same face-map idiom used in `remove_isolated_components!`). Using `face_key` keeps
this allocation-free, which matters because it runs once per solid tet in
`resolve_surface_candidates!`.
"""
function tetrahedron_faces(tet::Vector{Int64})
    a, b, c, d = tet[1], tet[2], tet[3], tet[4]
    return (
        face_key(a, b, c),
        face_key(a, b, d),
        face_key(a, c, d),
        face_key(b, c, d),
    )
end

"""
    resolve_surface_candidates!(mesh, solid_tets, candidates, scheme) -> Vector{Vector{Int64}}

Decide which deferred "quadruple-zero" tetrahedra (all four vertices on the surface) to keep,
following Labelle's surface-fidelity heuristic (isosurface-stuffing paper, §3.4):

  1. Discard a candidate that is inverted, or whose dihedral angles fall outside the [min, max]
     bounds -- such tetrahedra are too flat to improve surface fidelity.
  2. Of the well-shaped survivors, count how many of the four faces adjoin the solid mesh:
       * all four faces adjoin  -> RETAIN  (the tet fills a tetrahedral pocket in the boundary);
       * no face adjoins        -> DISCARD (an isolated "bubble");
       * otherwise              -> decide by the SDF sign at the centroid (inside -> keep).

This is a principled superset of quartet's `remove_exterior_tets` (which uses only the centroid
test) and is the PRIMARY mechanism for removing bubbles; `remove_isolated_components!` is left
afterwards only as a safety net that should now find almost nothing.

# Arguments
- `mesh::BlockMesh`: provides vertex coordinates and the SDF (read-only here)
- `solid_tets::Vector{Vector{Int64}}`: the already-decided interior mesh (no quadruple-zero tets)
- `candidates::Vector{Vector{Int64}}`: the deferred surface-only tetrahedra
- `scheme::String`: discretization scheme, used to fetch the dihedral-angle bounds

# Returns
- The subset of `candidates` to append to the mesh.
"""
function resolve_surface_candidates!(
    mesh::BlockMesh,
    solid_tets::Vector{Vector{Int64}},
    candidates::Vector{Vector{Int64}},
    scheme::String,
)::Vector{Vector{Int64}}
    retained = Vector{Vector{Int64}}()
    isempty(candidates) && return retained

    # Dihedral-angle bounds (Labelle §3.4): a quadruple-zero tet has all four vertices warped, so
    # reject it if any interior dihedral falls outside the accepted range (currently 10 .. 140 deg).
    bounds = create_warping_params(scheme)
    min_bound = bounds.min_dihedral_angle
    max_bound = bounds.max_dihedral_angle

    # (1) Apply the cheap, adjacency-INDEPENDENT filter first: drop inverted or too-flat
    # candidates, exactly as before. The raw A15 lattice tet is positively oriented and the
    # per-cell map preserves orientation, so a warped tet with non-positive signed volume has
    # been turned inside-out -> discard it (do NOT flip it); compute_dihedral_angle_range then
    # rejects the too-flat survivors. Only the survivors reach the adjacency test below, and we
    # record their four face keys while we have them.
    survivors = Vector{Vector{Int64}}()
    survivor_faces = Vector{NTuple{4,NTuple{3,Int}}}()
    for tet in candidates
        if !check_tetrahedron_orientation(mesh, tet)
            continue
        end
        min_angle, max_angle = compute_dihedral_angle_range(mesh, tet)
        if min_angle < min_bound || max_angle > max_bound
            continue
        end
        push!(survivors, tet)
        push!(survivor_faces, tetrahedron_faces(tet))
    end
    isempty(survivors) && return retained

    # (2) Count how many solid tets carry each survivor face. Seed the map with ONLY the
    # survivors' faces (value 0), then make a SINGLE pass over the solid mesh incrementing just
    # those faces. A candidate face "adjoins" the solid mesh exactly when some solid tet carries
    # it (value >= 1) -- the same presence test the old code did with haskey() over a full-mesh
    # map, except this Dict holds <= 4 * n_survivors entries instead of one per face of all ~2M
    # solid tets (building that full map is what dominated the slice cost). Candidates are not in
    # solid_tets, so a face shared by two candidates stays at 0 unless a solid tet also has it --
    # matching the old behaviour exactly.
    solid_face_count = Dict{NTuple{3,Int},Int}()
    for faces in survivor_faces
        for face in faces
            solid_face_count[face] = 0
        end
    end
    for tet in solid_tets
        for face in tetrahedron_faces(tet)
            if haskey(solid_face_count, face)
                solid_face_count[face] += 1
            end
        end
    end

    # (3) Decide each survivor with the SAME rule as before: count how many of its four faces
    # adjoin the solid mesh, then 4 -> retain, 0 -> discard (bubble), else centroid SDF sign.
    for k in eachindex(survivors)
        tet = survivors[k]
        adjoining = 0
        for face in survivor_faces[k]
            if solid_face_count[face] >= 1
                adjoining += 1
            end
        end

        if adjoining == 4
            push!(retained, tet)            # fills a tetrahedral pocket -> keep
        elseif adjoining == 0
            continue                        # isolated bubble -> discard
        else
            # Ambiguous: keep only if the centroid lies inside the geometry (quartet's test).
            centroid =
                (mesh.X[tet[1]] + mesh.X[tet[2]] + mesh.X[tet[3]] + mesh.X[tet[4]]) / 4.0
            if eval_sdf(mesh, centroid) < 0.0
                push!(retained, tet)
            end
        end
    end

    return retained
end

const RED = "\e[31m"
const BOLD = "\e[1m"
const RESET = "\e[0m"

"""
    remove_inverted_elements!(mesh::BlockMesh)

Improves mesh quality by fixing inverted tetrahedra and removing degenerate elements.

This function:
1. Fixes elements with negative orientation (inverted elements) by swapping two vertices.
   The orientation SIGN is decided EXACTLY by `is_positively_oriented` (no tolerance).
2. Removes elements with near-zero volume (degenerate elements). This is the only place a
   tolerance is kept, and it gates the volume MAGNITUDE, not the sign.
3. Leaves only the valid elements in `mesh.IEN`. It does NOT rebuild the inverse connectivity
   (INE) or compact orphaned nodes -- the caller refreshes connectivity right after (the pipeline
   and the test via `update_connectivity!`, `warp_mesh_by_planes_sdf!` via its own `create_INE!`),
   so rebuilding INE here would only be thrown away.

Returns the modified mesh.
"""
function remove_inverted_elements!(mesh::BlockMesh)
    @info "Fixing elements orientation..."
    # Tolerance used ONLY to drop near-zero-volume (degenerate) elements. The orientation SIGN
    # itself is decided exactly by is_positively_oriented, with no tolerance.
    volume_tolerance = mesh.grid_tol * 1e-6

    # Track statistics for reporting
    fixed_elements = 0
    failed_fixes = 0
    zero_volume_elements = 0

    # Process elements: fix inverted elements, remove near-zero volume elements
    valid_elements = Vector{Vector{Int}}()
    sizehint!(valid_elements, length(mesh.IEN))

    for tet in mesh.IEN
        # Skip degenerate elements with duplicate vertices
        if !four_distinct(tet[1], tet[2], tet[3], tet[4])
            zero_volume_elements += 1
            continue
        end

        # Float signed volume, used here ONLY for the near-zero-volume MAGNITUDE test.
        a = mesh.X[tet[2]] - mesh.X[tet[1]]
        b = mesh.X[tet[3]] - mesh.X[tet[1]]
        c = mesh.X[tet[4]] - mesh.X[tet[1]]
        if abs(dot(a, cross(b, c))) <= volume_tolerance
            zero_volume_elements += 1
            continue
        end

        # Orientation SIGN decided exactly. A non-degenerate tetrahedron is either already
        # positively oriented, or a single vertex swap (3 <-> 4) flips it to positive: one
        # transposition exactly negates a non-zero determinant, so the old fallback strategies
        # (which only existed to paper over floating-point noise) are no longer needed.
        if is_positively_oriented(mesh, tet)
            push!(valid_elements, tet)
        else
            tet_fixed = copy(tet)
            tet_fixed[3], tet_fixed[4] = tet_fixed[4], tet_fixed[3]
            if is_positively_oriented(mesh, tet_fixed)
                push!(valid_elements, tet_fixed)
                fixed_elements += 1
            else
                # Unreachable for a tet that passed the volume test above (it cannot be exactly
                # coplanar); kept as a safety net.
                failed_fixes += 1
            end
        end
    end

    # Keep only the valid elements. INE is intentionally NOT rebuilt here (see the docstring):
    # every caller refreshes connectivity right after, so a rebuild now would be discarded.
    mesh.IEN = valid_elements

    # Report statistics
    println("  Fixed orientation of $fixed_elements inverted elements")
    if failed_fixes != 0
        println(
            "  Failed to fix orientation of $(RED)$(BOLD)$(failed_fixes)$(RESET) elements",
        )
    end
    println("  Removing $zero_volume_elements elements with near-zero volume")

    return mesh
end
