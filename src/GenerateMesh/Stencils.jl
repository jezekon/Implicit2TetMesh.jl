# Functions for tetrahedral mesh slicing along an isosurface

"""
    slice_ambiguous_tetrahedra!(mesh::BlockMesh, scheme::String)

Slice tetrahedra that stick out across the isosurface (SDF zero level set), trimming each one
back to the interior region (SDF <= 0). This is a direct port of quartet's `trim_spikes`
(make_tet_mesh.cpp:221-346). Every spike tetrahedron is replaced by smaller tetrahedra whose
vertices lie either on the original lattice or on the linear edge cut points produced by
`cut_edge!`.

The "quadruple-zero" tetrahedra (all four vertices on the surface) are handled separately in a
SECOND pass by `resolve_surface_candidates!`, following Labelle's surface-fidelity heuristic
(§3.4). They have to be deferred because their fate depends on whether their faces adjoin the rest
of the solid mesh, which is only known once the whole solid mesh has been built.

# Arguments
- `mesh::BlockMesh`: The mesh to process
- `scheme::String`: Discretization scheme (only "A15" is supported)
"""
function slice_ambiguous_tetrahedra!(mesh::BlockMesh, scheme::String)
    @info "Slicing tetrahedra using trim_spikes logic..."

    cut_map = Dict{Tuple{Int,Int},Int}()
    solid_tets = Vector{Vector{Int64}}()          # definitely-interior output tets (no quadruple-zero)
    surface_candidates = Vector{Vector{Int64}}()  # all-on-surface tets, resolved in the second pass
    sizehint!(solid_tets, length(mesh.IEN))

    current_IEN = mesh.IEN
    mesh.IEN = Vector{Vector{Int64}}()

    # First pass: trim every tetrahedron. Crossing and fully-interior tets are emitted straight
    # into the solid mesh; all-on-surface ("quadruple-zero") tets are collected as candidates.
    for tet in current_IEN
        resulting_tets = apply_stencil_trim_spikes!(mesh, tet, cut_map, surface_candidates)
        for nt in resulting_tets
            push!(solid_tets, nt)
        end
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
    cut_edge!(i, j, mesh, node_sdf, cut_map) -> Int

Find or create the vertex where the isosurface crosses the edge `(i, j)`, placed at the LINEAR
interpolation point

    X_cut = (1 - alpha) * X_i + alpha * X_j,    alpha = phi_i / (phi_i - phi_j)  in (0, 1),

exactly as in quartet's `cut_edge` (make_tet_mesh.cpp:196-217). Using the linear point (instead
of projecting onto the true isosurface) keeps the trimmed tetrahedra well shaped and is
consistent with the edge-based `warp!`. Each edge is cut only once: the new vertex index is
cached in `cut_map` (keyed by the canonical edge), so all tetrahedra sharing the edge reuse the
same vertex and the boundary stays watertight.

If an endpoint already lies on the surface (`|phi| < tol`) that endpoint index is returned
instead of creating a new vertex.
"""
function cut_edge!(
    i::Int,
    j::Int,
    mesh::BlockMesh,
    node_sdf::Vector{Float64},
    cut_map::Dict{Tuple{Int,Int},Int},
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

    # Linear cut point along the edge (quartet cut_edge): alpha is the fraction from i to the
    # zero crossing, and the point is the straight-line interpolation between the endpoints.
    alpha = sdf_i / (sdf_i - sdf_j)
    cut_point = (1.0 - alpha) * mesh.X[i] + alpha * mesh.X[j]

    push!(mesh.X, cut_point)
    push!(mesh.node_sdf, 0.0)   # the cut point lies exactly on the surface
    new_index = length(mesh.X)
    cut_map[edge] = new_index
    return new_index
end

"""
    check_tetrahedron_orientation(mesh::BlockMesh, tet::Vector{Int})

Check if a tetrahedron has positive orientation (positive Jacobian determinant).
Returns true for correctly oriented tetrahedra, false for inverted ones.
"""
function check_tetrahedron_orientation(mesh::BlockMesh, tet::Vector{Int})
    # Validate indices
    if length(tet) != 4 || any(i -> i <= 0 || i > length(mesh.X), tet)
        @warn "Invalid tetrahedron indices for orientation check: $tet"
        return false
    end

    # Get tetrahedron vertices
    vertices = [mesh.X[tet[i]] for i = 1:4]

    # Calculate edge vectors from first vertex
    a = vertices[2] - vertices[1]
    b = vertices[3] - vertices[1]
    c = vertices[4] - vertices[1]

    # Calculate Jacobian determinant (proportional to signed volume)
    det_value = dot(a, cross(b, c))

    # Positive determinant indicates correct orientation
    return det_value > 1e-12
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
    apply_stencil_trim_spikes!(mesh, tet, cut_map, surface_candidates)

Apply the trimming algorithm to a single tetrahedron that crosses the isosurface.
Returns new tetrahedra that accurately represent the interior region (SDF <= 0).
The SDF sign pattern is classified case-by-case, ported 1:1 from quartet `trim_spikes`
(make_tet_mesh.cpp:221-346).

A "quadruple-zero" tetrahedron (all four vertices on the surface) is NOT decided here: it is
appended to `surface_candidates` and resolved later by `resolve_surface_candidates!` (Labelle
§3.4), because its fate depends on the surrounding solid mesh.
"""
function apply_stencil_trim_spikes!(
    mesh::BlockMesh,
    tet::Vector{Int64},
    cut_map::Dict{Tuple{Int,Int},Int},
    surface_candidates::Vector{Vector{Int64}},
)::Vector{Vector{Int64}}

    # Get SDF values for tetrahedron nodes
    node_indices = tet
    node_sdf = [mesh.node_sdf[idx] for idx in node_indices]
    tol = mesh.grid_tol

    # --- Special cases handling ---
    # Case 1: All nodes outside (SDF > tol) - discard tetrahedron (all-outside case)
    if all(s -> s > tol, node_sdf)
        return Vector{Vector{Int64}}()
    end

    # Case 2: All nodes inside (SDF < -tol) - keep original (all-inside case)
    if all(s -> s < -tol, node_sdf)
        return [tet]
    end

    # Case 3: All nodes on the surface (|SDF| <= tol) - "quadruple-zero" tetrahedron. Its fate is
    # ambiguous (Labelle §3.4) and depends on the surrounding solid mesh, so defer it: collect it
    # as a candidate and emit nothing now. resolve_surface_candidates! decides it afterwards.
    if all(s -> abs(s) <= tol, node_sdf)
        push!(surface_candidates, tet)
        return Vector{Vector{Int64}}()
    end

    # --- General case: Tetrahedron crossing the isosurface ---
    # Sort vertices by SDF value (largest first, i.e. most-outside vertex first)
    vert_data = [(node_sdf[i], node_indices[i]) for i = 1:4]
    p = [1, 2, 3, 4]
    flipped = false

    # Comparison function for consistent vertex ordering, matching quartet trim_spikes.
    # Descending by SDF; ties are broken so the LARGER original node index is treated as the
    # more-outside vertex (the same direction quartet uses when it swaps on `index <`). Because
    # this is a global total order on (SDF, index), face-adjacent tetrahedra sort a shared quad
    # identically and therefore split it along the same diagonal -- this is what keeps the
    # boundary watertight.
    less_than(idx1, idx2) =
        (vert_data[idx1][1] > vert_data[idx2][1]) || (
            vert_data[idx1][1] == vert_data[idx2][1] &&
            vert_data[idx1][2] > vert_data[idx2][2]
        )

    # Sort vertices while tracking orientation flips
    if less_than(p[2], p[1])
        p[1], p[2] = p[2], p[1];
        flipped = !flipped
    end
    if less_than(p[4], p[3])
        p[3], p[4] = p[4], p[3];
        flipped = !flipped
    end
    if less_than(p[3], p[1])
        p[1], p[3] = p[3], p[1];
        flipped = !flipped
    end
    if less_than(p[4], p[2])
        p[2], p[4] = p[4], p[2];
        flipped = !flipped
    end
    if less_than(p[3], p[2])
        p[2], p[3] = p[3], p[2];
        flipped = !flipped
    end

    # Sorted vertices (s has highest SDF = most outside, p has lowest = most inside)
    s_idx, r_idx, q_idx, p_idx = [vert_data[pi][2] for pi in p]
    sdf_s, sdf_r, sdf_q, sdf_p = [vert_data[pi][1] for pi in p]

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

    new_tets = Vector{Vector{Int64}}()

    # Helper function to add tetrahedron with orientation check
    function add_tet!(t::Vector{Int})
        # Skip degenerate cases (duplicate vertices)
        if length(Set(t)) != 4
            return
        end

        # Fix orientation and add to result
        fix_tetrahedron_orientation!(mesh, t)
        if !check_tetrahedron_orientation(mesh, t)
            @warn "Tetrahedron $t still has incorrect orientation after fix. Skipping."
            return
        end
        push!(new_tets, t)
    end

    # --- Case analysis based on SDF sign patterns --- (S >= R >= Q >= P)

    # Entirely outside or on the surface (quartet's vphi[s]==0 branch: +++0 / ++00 / +000).
    # The most-inside vertex is already on the surface, so no interior region remains -> discard.
    if is_p_surface
        return Vector{Vector{Int64}}()

        # Surface/interior tetrahedron with no outside vertex: ZZZP, ZZPP, ZPPP. quartet leaves
        # these untouched (all phi <= 0), so we keep them as-is.
    elseif is_s_surface && !is_r_outside && !is_q_outside && is_p_inside
        return [tet]

        # Case NNNP (quartet +++-): three nodes outside, one inside
    elseif is_s_outside && is_r_outside && is_q_outside && is_p_inside
        # Cut three edges from outside nodes to inside node
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        qp = cut_edge!(q_idx, p_idx, mesh, mesh.node_sdf, cut_map)

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
            sr = cut_edge!(s_idx, r_idx, mesh, mesh.node_sdf, cut_map)
            sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
            sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)

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
            return Vector{Vector{Int64}}()
        end

        # Case NNPP (quartet ++--): two nodes outside, two inside
    elseif is_r_outside && is_q_inside
        # Cut all four edges crossing the isosurface
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rq = cut_edge!(r_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)

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
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)

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
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)

        # Create one tetrahedron
        if flipped
            add_tet!([p_idx, q_idx, rp, sp])
        else
            add_tet!([q_idx, p_idx, rp, sp])
        end

        # Case NZZP (quartet +00-): one outside, two on surface, one inside
    elseif is_s_outside && is_r_surface && is_q_surface && is_p_inside
        # Cut edge from outside node to inside node
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)

        # Create one tetrahedron
        if flipped
            add_tet!([q_idx, r_idx, p_idx, sp])
        else
            add_tet!([r_idx, q_idx, p_idx, sp])
        end

    else
        @warn "Unexpected SDF pattern in apply_stencil_trim_spikes!: s=$sdf_s, r=$sdf_r, q=$sdf_q, p=$sdf_p. Indices: s=$s_idx, r=$r_idx, q=$q_idx, p=$p_idx. Tet: $tet. Flipped: $flipped"
        return Vector{Vector{Int64}}() # Discard in case of unexpected pattern
    end

    return new_tets
end

"""
    tetrahedron_faces(tet) -> NTuple{4, NTuple{3,Int}}

Return the four triangular faces of a tetrahedron as sorted node-index triples. Sorting makes the
key orientation-independent, so a face shared by two tetrahedra maps to the same key (the same
face-map idiom used in `remove_isolated_components!`).
"""
function tetrahedron_faces(tet::Vector{Int64})
    a, b, c, d = tet[1], tet[2], tet[3], tet[4]
    return (
        Tuple(sort([a, b, c])),
        Tuple(sort([a, b, d])),
        Tuple(sort([a, c, d])),
        Tuple(sort([b, c, d])),
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

    # Build the face -> incidence-count map of the solid mesh. A candidate face "adjoins" the
    # solid mesh exactly when it is present here (reuse of the remove_isolated_components! idiom).
    solid_face_count = Dict{NTuple{3,Int},Int}()
    for tet in solid_tets
        for face in tetrahedron_faces(tet)
            solid_face_count[face] = get(solid_face_count, face, 0) + 1
        end
    end

    # Dihedral-angle bounds. A quadruple-zero tet has all four vertices warped, so we use the
    # stricter NNZZ bounds from create_warping_params (currently 10 deg .. 140 deg).
    params = create_warping_params(scheme, mesh.grid_step)
    min_bound = params.nnzz.min_dihedral_angle
    max_bound = params.nnzz.max_dihedral_angle

    for tet in candidates
        # (1) Reject inverted or badly-shaped tets. The raw A15 lattice tet is positively
        # oriented and the per-cell map preserves orientation, so a warped tet with non-positive
        # signed volume has been turned inside-out -> discard it (do NOT flip it).
        # compute_dihedral_angle_range then rejects the too-flat survivors.
        if !check_tetrahedron_orientation(mesh, tet)
            continue
        end
        min_angle, max_angle = compute_dihedral_angle_range(mesh, tet)
        if min_angle < min_bound || max_angle > max_bound
            continue
        end

        # (2) Count how many of the four faces adjoin the solid mesh.
        adjoining = 0
        for face in tetrahedron_faces(tet)
            if haskey(solid_face_count, face)
                adjoining += 1
            end
        end

        if adjoining == 4
            push!(retained, tet)            # fills a tetrahedral pocket -> keep
        elseif adjoining == 0
            continue                        # isolated bubble -> discard
        else
            # (3) Ambiguous: keep only if the centroid lies inside the geometry (quartet's test).
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
1. Attempts to fix elements with negative Jacobian determinant (inverted elements)
   by reordering their vertices to achieve positive orientation
2. Removes elements with near-zero volume (degenerate elements)
3. Updates mesh connectivity to remove orphaned nodes
4. Rebuilds the inverse node-to-element connectivity (INE)

Returns the modified mesh.
"""
function remove_inverted_elements!(mesh::BlockMesh)
    @info "Fixing elements orientation..."
    # Set tolerance for identifying near-zero volumes
    volume_tolerance = mesh.grid_tol * 1e-6

    # Track statistics for reporting
    fixed_elements = 0
    failed_fixes = 0
    zero_volume_elements = 0

    # Process elements: fix inverted elements, remove near-zero volume elements
    valid_elements = Vector{Vector{Int}}()
    sizehint!(valid_elements, length(mesh.IEN))

    for (elem_idx, tet) in enumerate(mesh.IEN)
        # Skip degenerate elements with duplicate vertices
        if length(Set(tet)) != 4
            zero_volume_elements += 1
            continue
        end

        # Calculate determinant (proportional to signed volume)
        vertices = SVector{4,SVector{3,Float64}}(
            mesh.X[tet[1]],
            mesh.X[tet[2]],
            mesh.X[tet[3]],
            mesh.X[tet[4]],
        )
        a = vertices[2] - vertices[1]
        b = vertices[3] - vertices[1]
        c = vertices[4] - vertices[1]
        det_value = dot(a, cross(b, c))

        # Remove elements with near-zero volume
        if abs(det_value) <= volume_tolerance
            zero_volume_elements += 1
            continue
        end

        # Fix elements with negative determinant
        if det_value < 0
            # Strategy 1: Swap vertices 3 and 4
            tet_copy = copy(tet)
            tet_copy[3], tet_copy[4] = tet_copy[4], tet_copy[3]

            # Check if fix worked
            new_vertices = SVector{4,SVector{3,Float64}}(
                mesh.X[tet_copy[1]],
                mesh.X[tet_copy[2]],
                mesh.X[tet_copy[3]],
                mesh.X[tet_copy[4]],
            )
            new_a = new_vertices[2] - new_vertices[1]
            new_b = new_vertices[3] - new_vertices[1]
            new_c = new_vertices[4] - new_vertices[1]
            new_det = dot(new_a, cross(new_b, new_c))

            if new_det > volume_tolerance
                # Fix successful - element has positive volume and is not near-zero
                push!(valid_elements, tet_copy)
                fixed_elements += 1
                continue
            end

            # Strategy 2: Swap vertices 1 and 2
            tet_copy = copy(tet)
            tet_copy[1], tet_copy[2] = tet_copy[2], tet_copy[1]

            new_vertices = SVector{4,SVector{3,Float64}}(
                mesh.X[tet_copy[1]],
                mesh.X[tet_copy[2]],
                mesh.X[tet_copy[3]],
                mesh.X[tet_copy[4]],
            )
            new_a = new_vertices[2] - new_vertices[1]
            new_b = new_vertices[3] - new_vertices[1]
            new_c = new_vertices[4] - new_vertices[1]
            new_det = dot(new_a, cross(new_b, new_c))

            if new_det > volume_tolerance
                push!(valid_elements, tet_copy)
                fixed_elements += 1
                continue
            end

            # Strategy 3: Swap vertices 2 and 3
            tet_copy = copy(tet)
            tet_copy[2], tet_copy[3] = tet_copy[3], tet_copy[2]

            new_vertices = SVector{4,SVector{3,Float64}}(
                mesh.X[tet_copy[1]],
                mesh.X[tet_copy[2]],
                mesh.X[tet_copy[3]],
                mesh.X[tet_copy[4]],
            )
            new_a = new_vertices[2] - new_vertices[1]
            new_b = new_vertices[3] - new_vertices[1]
            new_c = new_vertices[4] - new_vertices[1]
            new_det = dot(new_a, cross(new_b, new_c))

            if new_det > volume_tolerance
                push!(valid_elements, tet_copy)
                fixed_elements += 1
                continue
            end

            # All fix attempts failed
            failed_fixes += 1
        else
            # Element already has positive determinant
            push!(valid_elements, tet)
        end
    end

    # Update connectivity
    mesh.IEN = valid_elements
    create_INE!(mesh)                  # Creates inverse connectivity (mesh.INE)

    # Report statistics before connectivity update
    println("  Fixed orientation of $fixed_elements inverted elements")
    if failed_fixes != 0
        println(
            "  Failed to fix orientation of $(RED)$(BOLD)$(failed_fixes)$(RESET) elements",
        )
    end
    println("  Removing $zero_volume_elements elements with near-zero volume")

    return mesh
end
