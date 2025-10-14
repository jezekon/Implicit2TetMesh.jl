# Functions for tetrahedral mesh slicing along an isosurface
# Implements the "trim spikes" algorithm for accurate surface representation

"""
    slice_ambiguous_tetrahedra!(mesh::BlockMesh)

Slice tetrahedra crossing the isosurface (SDF zero level set).
Identifies elements crossing the boundary and replaces them with smaller
tetrahedra that accurately represent the surface.
"""
function slice_ambiguous_tetrahedra!(mesh::BlockMesh)
    @info "Slicing tetrahedra using trim_spikes logic (with element splitting)..."
    new_IEN = Vector{Vector{Int64}}()
    sizehint!(new_IEN, length(mesh.IEN)) # Pre-allocate for performance

    # Cache for storing cut edge results to avoid duplicate cutting
    cut_map = Dict{Tuple{Int,Int},Int}()

    original_node_count = length(mesh.X)
    original_tet_count = length(mesh.IEN)

    # Process original connectivity while building new connectivity
    current_IEN = mesh.IEN
    mesh.IEN = Vector{Vector{Int64}}() # Clear current IEN temporarily

    for tet in current_IEN
        # Apply the trim_spikes stencil to each tetrahedron
        resulting_tets = apply_stencil_trim_spikes!(mesh, tet, cut_map)

        # Add valid resulting tetrahedra to the new connectivity list
        for nt in resulting_tets
            push!(new_IEN, nt)
        end
    end

    mesh.IEN = new_IEN # Update mesh connectivity
    new_tet_count = length(mesh.IEN)

    println(
        "  After trim_spikes slicing: $(new_tet_count) tetrahedra (removed $(original_tet_count - new_tet_count))",
    )
end

"""
    cut_edge!(i::Int, j::Int, mesh::BlockMesh, node_sdf::Vector{Float64}, cut_map::Dict{Tuple{Int, Int}, Int})::Int

Calculate the intersection point between an edge (i,j) and the isosurface.
Returns the node index at the intersection point, creating a new node if needed.
Handles special cases where nodes are already on the surface.
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

    # Handle cases where nodes are already on the surface
    if is_i_zero && is_j_zero
        return min(i, j)
    elseif is_i_zero
        return i
    elseif is_j_zero
        return j
    end

    # Ensure nodes are on opposite sides
    if sign(sdf_i) == sign(sdf_j)
        error(
            "cut_edge! called with nodes on same side: i=$i (sdf=$sdf_i), j=$j (sdf=$sdf_j)",
        )
    end

    # Canonical edge representation
    edge = (min(i, j), max(i, j))

    if haskey(cut_map, edge)
        return cut_map[edge]
    end

    # Get node positions
    pos_i = mesh.X[i]
    pos_j = mesh.X[j]

    new_index = interpolate_zero(pos_i, pos_j, sdf_i, sdf_j, mesh)

    # Cache result
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
    warp_node_to_surface(mesh::BlockMesh, position::SVector{3,Float64}, current_sdf::Float64)

Warp a node with negative SDF to the zero isosurface using gradient descent.
Returns the new position on the surface.
"""
function warp_node_to_surface(
    mesh::BlockMesh,
    position::SVector{3,Float64},
    current_sdf::Float64,
    max_iter::Int = 20,
)::SVector{3,Float64}

    tol = mesh.grid_tol
    current_pos = position

    # Newton iteration to find zero level set
    for iter = 1:max_iter
        f = eval_sdf(mesh, current_pos)

        # Close enough to surface
        if abs(f) < tol
            break
        end

        # Compute gradient
        grad = compute_gradient(mesh, current_pos)
        norm_grad_sq = sum(abs2, grad)

        # Avoid division by zero
        if norm_grad_sq < 1e-16
            @warn "Gradient too small during warp, stopping at SDF = $f"
            break
        end

        # Newton step: move along gradient toward zero level
        dp = (f / norm_grad_sq) * grad
        current_pos -= dp
    end

    # Verify we reached the surface
    final_sdf = eval_sdf(mesh, current_pos)
    if abs(final_sdf) > tol * 100
        @warn "Warp did not converge to surface: final SDF = $final_sdf"
    end

    return current_pos
end

"""
    add_warped_node!(mesh::BlockMesh, position::SVector{3,Float64}, cut_map::Dict)

Add a new node to the mesh at the given position with SDF = 0.
Uses cut_map to avoid creating duplicate nodes at the same location.
"""
function add_warped_node!(
    mesh::BlockMesh,
    position::SVector{3,Float64},
    cut_map::Dict{Tuple{Int,Int},Int},
)::Int

    tol = mesh.grid_tol

    # Quantize position to handle numerical precision
    p_key = quantize(position, tol)

    # Check if node already exists at this location
    if haskey(mesh.node_hash, p_key)
        return mesh.node_hash[p_key]
    end

    # Create new node
    push!(mesh.X, position)
    push!(mesh.node_sdf, 0.0)  # Surface node
    new_index = length(mesh.X)
    mesh.node_hash[p_key] = new_index

    return new_index
end

"""
    apply_stencil_trim_spikes!(mesh::BlockMesh, tet::Vector{Int64}, cut_map::Dict{Tuple{Int, Int}, Int})

Apply trimming algorithm to a tetrahedron that crosses the isosurface.
Returns new tetrahedra that accurately represent the interior region (SDF ≥ 0).
Uses case-by-case analysis based on SDF sign patterns for robust slicing.
"""
function apply_stencil_trim_spikes!(
    mesh::BlockMesh,
    tet::Vector{Int64},
    cut_map::Dict{Tuple{Int,Int},Int},
)::Vector{Vector{Int64}}

    # Get SDF values for tetrahedron nodes
    node_indices = tet
    node_sdf = [mesh.node_sdf[idx] for idx in node_indices]
    tol = mesh.grid_tol

    # --- Special cases handling ---
    # Case 1: All nodes outside (SDF < -tol) - discard tetrahedron (NNNN case)
    if all(s -> s < -tol, node_sdf)
        return Vector{Vector{Int64}}()
    end

    # Case 2: All nodes inside or on surface (SDF ≥ tol) - keep original (PPPP case)
    if all(s -> s > tol, node_sdf)
        return [tet]
    end

    # --- General case: Tetrahedron crossing the isosurface ---
    # Sort vertices by SDF value (smallest first)
    vert_data = [(node_sdf[i], node_indices[i]) for i = 1:4]
    p = [1, 2, 3, 4]
    flipped = false

    # Comparison function for consistent vertex ordering
    less_than(idx1, idx2) =
        (vert_data[idx1][1] < vert_data[idx2][1]) || (
            vert_data[idx1][1] == vert_data[idx2][1] &&
            vert_data[idx1][2] < vert_data[idx2][2]
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

    # Sorted vertices (s has lowest SDF, p has highest)
    s_idx, r_idx, q_idx, p_idx = [vert_data[pi][2] for pi in p]
    sdf_s, sdf_r, sdf_q, sdf_p = [vert_data[pi][1] for pi in p]

    # Check if SDF values are properly sorted
    if !(sdf_s <= sdf_r <= sdf_q <= sdf_p)
        @warn "SDF values not sorted: s=$sdf_s, r=$sdf_r, q=$sdf_q, p=$sdf_p"
    end

    # Case 3: All nodes on surface (abs(SDF) < tol) - keep original (ZZZZ case)
    if all(s -> abs(s) <= tol, node_sdf)
        centroid = (mesh.X[s_idx] + mesh.X[r_idx] + mesh.X[q_idx] + mesh.X[p_idx]) / 4.0
        # println("ZZZZ case, nodes sdf: $(node_sdf)")
        if eval_sdf(mesh, centroid) > 0.0
            # This tetrahedron represents the surface, keep it
            return [tet]
        else
            # Discard if centroid is outside
            return Vector{Vector{Int64}}()
        end
    end

    # Classify vertices by SDF sign (N=negative, P=positive, Z=zero)
    is_s_neg = sdf_s < -tol
    is_r_neg = sdf_r < -tol
    is_q_neg = sdf_q < -tol
    is_p_neg = sdf_p < -tol

    # is_s_pos = sdf_s > tol
    is_p_pos = sdf_p > tol
    is_q_pos = sdf_q > tol
    is_r_pos = sdf_r > tol

    is_s_zero = abs(sdf_s) <= tol
    is_r_zero = abs(sdf_r) <= tol
    is_q_zero = abs(sdf_q) <= tol
    is_p_zero = abs(sdf_p) <= tol

    new_tets = Vector{Vector{Int64}}()

    # Helper function to add tetrahedron with orientation check
    function add_tet!(t::Vector{Int})
        # Skip degenerate cases (duplicate vertices)
        if length(Set(t)) != 4
            # @warn "Degenerate tetrahedron generated (duplicate nodes): $t. Skipping."
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

    # --- Case analysis based on SDF sign patterns --- (S=<R=<Q=<P)
    if is_s_neg && is_p_zero # NNNZ, NNZZ, NZZZ
        # Compute centroid to check if element should be preserved
        centroid = (mesh.X[s_idx] + mesh.X[r_idx] + mesh.X[q_idx] + mesh.X[p_idx]) / 4.0
        centroid_sdf = eval_sdf(mesh, centroid)

        # If centroid is inside (positive SDF), we need to handle this carefully
        if centroid_sdf > tol
            # Identify which specific case we're in
            if is_s_neg && is_r_neg && is_q_neg && is_p_zero
                # NNNZ - three nodes outside, one on surface
                # This is very thin, likely discard
                return Vector{Vector{Int64}}()

            elseif is_s_neg && is_r_neg && is_q_zero && is_p_zero
                # NNZZ - two nodes outside, two on surface
                # Warp negative nodes to surface and create new element

                # Warp both negative nodes to surface
                s_warped = warp_node_to_surface(mesh, mesh.X[s_idx], sdf_s)
                r_warped = warp_node_to_surface(mesh, mesh.X[r_idx], sdf_r)

                # Add new warped nodes to mesh
                s_new = add_warped_node!(mesh, s_warped, cut_map)
                r_new = add_warped_node!(mesh, r_warped, cut_map)

                # Create new tetrahedron with warped nodes
                new_tet = [s_new, r_new, q_idx, p_idx]

                # Fix orientation if needed
                fix_tetrahedron_orientation!(mesh, new_tet)

                # Verify orientation
                if check_tetrahedron_orientation(mesh, new_tet)
                    return [new_tet]
                else
                    @warn "NNZZ warped tetrahedron has bad orientation, discarding"
                    return Vector{Vector{Int64}}()
                end

            elseif is_s_neg && is_r_zero && is_q_zero && is_p_zero
                # NZZZ - one node outside, three on surface
                # Warp the single negative node to surface

                # Warp negative node to surface
                s_warped = warp_node_to_surface(mesh, mesh.X[s_idx], sdf_s)

                # Add new warped node to mesh
                s_new = add_warped_node!(mesh, s_warped, cut_map)

                # Create new tetrahedron with warped node
                new_tet = [s_new, r_idx, q_idx, p_idx]

                # Fix orientation if needed
                fix_tetrahedron_orientation!(mesh, new_tet)

                # Verify orientation
                if check_tetrahedron_orientation(mesh, new_tet)
                    return [new_tet]
                else
                    @warn "NZZZ warped tetrahedron has bad orientation, discarding"
                    return Vector{Vector{Int64}}()
                end
            end
        end

        # If centroid is negative or we didn't match specific cases, discard
        return Vector{Vector{Int64}}()

        # Case (Z,!N,!N,P): Surface node with all others on surface or inside
    elseif is_s_zero && !is_r_neg && !is_q_neg && is_p_pos # ZZZP, ZZPP, ZPPP
        # println("ZZZP, ZZPP, ZPPP")
        return [tet]

        # Case NNNP: Three nodes outside, one inside
    elseif is_s_neg && is_r_neg && is_q_neg && is_p_pos
        # println("NNNP")
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

        # Case NPPP: One node outside, three inside
    elseif is_s_neg && is_r_pos
        # println("NPPP")
        if !is_r_neg && !is_q_neg # Confirm r, q, p are interior nodes
            # Cut edges from outside node to each interior node
            sr = cut_edge!(s_idx, r_idx, mesh, mesh.node_sdf, cut_map)
            sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
            sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)

            # Create three tetrahedra filling the interior region
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

        # Case NNPP: Two nodes outside, two inside
    elseif is_r_neg && is_q_pos
        # println("NNPP")
        # Cut all four edges crossing the isosurface
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rq = cut_edge!(r_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)

        # Create three tetrahedra for interior region
        if flipped
            add_tet!([p_idx, rq, q_idx, sq])
            add_tet!([p_idx, sp, sq, rp])
            add_tet!([p_idx, rp, sq, rq])
        else
            add_tet!([p_idx, q_idx, rq, sq])
            add_tet!([p_idx, sq, sp, rp])
            add_tet!([p_idx, sq, rp, rq])
        end

        # Case NZPP: One outside, one on surface, two inside
    elseif is_s_neg && is_r_zero && is_q_pos
        # println("NZPP")
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

        # Case NNZP: Two outside, one on surface, one inside
    elseif is_r_neg && is_q_zero && is_p_pos
        # println("NNZP")
        # Cut edges from outside nodes to inside node
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)

        # Create one tetrahedron
        if flipped
            add_tet!([p_idx, q_idx, rp, sp])
        else
            add_tet!([q_idx, p_idx, rp, sp])
        end

        # Case NZZP: One outside, two on surface, one inside
    elseif is_s_neg && is_r_zero && is_q_zero && is_p_pos
        # println("NZZP")
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
    fix_tetrahedra_orientation!(mesh::BlockMesh)

Checks and fixes the orientation of all tetrahedral elements in the mesh
by ensuring each element has a positive Jacobian determinant (positive volume).
Elements with negative determinants are fixed by reordering their nodes.

# Arguments
- `mesh::BlockMesh`: The mesh to check and fix
- `tol::Float64`: Tolerance for considering volume as zero (default: 1e-12)

# Returns
- Number of elements that were fixed
"""
function fix_tetrahedra_orientation!(mesh::BlockMesh, tol::Float64 = 1e-12)
    fixed_count = 0
    degenerate_count = 0
    failed_count = 0

    @info "Checking and fixing tetrahedra orientation..."

    for (elem_idx, tet) in enumerate(mesh.IEN)
        # Skip elements with duplicate nodes
        if length(Set(tet)) != 4
            degenerate_count += 1
            continue
        end

        # Get vertices of the tetrahedron
        vertices = [mesh.X[node_idx] for node_idx in tet]

        # Calculate edge vectors from first vertex
        a = vertices[2] - vertices[1]
        b = vertices[3] - vertices[1]
        c = vertices[4] - vertices[1]

        # Calculate Jacobian determinant (6 times the volume)
        det_value = dot(a, cross(b, c))

        # Check if element has negative or zero volume
        if det_value <= tol
            # Try to fix by swapping nodes 3 and 4
            mesh.IEN[elem_idx][3], mesh.IEN[elem_idx][4] =
                mesh.IEN[elem_idx][4], mesh.IEN[elem_idx][3]

            # Check if the fix worked
            vertices = [mesh.X[node_idx] for node_idx in mesh.IEN[elem_idx]]
            a = vertices[2] - vertices[1]
            b = vertices[3] - vertices[1]
            c = vertices[4] - vertices[1]
            new_det = dot(a, cross(b, c))

            if new_det > tol
                fixed_count += 1
            else
                # If that didn't work, try swapping nodes 1 and 2
                mesh.IEN[elem_idx][3], mesh.IEN[elem_idx][4] =
                    mesh.IEN[elem_idx][4], mesh.IEN[elem_idx][3] # Restore
                mesh.IEN[elem_idx][1], mesh.IEN[elem_idx][2] =
                    mesh.IEN[elem_idx][2], mesh.IEN[elem_idx][1]

                # Check again
                vertices = [mesh.X[node_idx] for node_idx in mesh.IEN[elem_idx]]
                a = vertices[2] - vertices[1]
                b = vertices[3] - vertices[1]
                c = vertices[4] - vertices[1]
                new_det = dot(a, cross(b, c))

                if new_det > tol
                    fixed_count += 1
                else
                    # Restore and try one more swap
                    mesh.IEN[elem_idx][1], mesh.IEN[elem_idx][2] =
                        mesh.IEN[elem_idx][2], mesh.IEN[elem_idx][1] # Restore
                    mesh.IEN[elem_idx][2], mesh.IEN[elem_idx][3] =
                        mesh.IEN[elem_idx][3], mesh.IEN[elem_idx][2]

                    vertices = [mesh.X[node_idx] for node_idx in mesh.IEN[elem_idx]]
                    a = vertices[2] - vertices[1]
                    b = vertices[3] - vertices[1]
                    c = vertices[4] - vertices[1]
                    new_det = dot(a, cross(b, c))

                    if new_det > tol
                        fixed_count += 1
                    else
                        failed_count += 1
                    end
                end
            end
        end
    end

    if degenerate_count > 0
        @warn "Found $degenerate_count elements with duplicate nodes"
    end

    if failed_count > 0
        @warn "Failed to fix orientation of $failed_count elements"
    end

    @info "Fixed orientation of $fixed_count tetrahedra"
    return fixed_count
end
