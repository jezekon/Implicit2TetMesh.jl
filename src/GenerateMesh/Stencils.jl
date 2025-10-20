# Functions for tetrahedral mesh slicing along an isosurface
# Implements the "trim spikes" algorithm for accurate surface representation

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
"""
struct CaseParams
    threshold_distance::Float64
    max_node_displacement::Float64
    min_volume_ratio::Float64
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
        nnzz = CaseParams(base_threshold, 0.2 * grid_step, 0.05)
        # NZZZ: Relaxed parameters (warping only 1 node)
        nzzz = CaseParams(2.0 * base_threshold, 0.4 * grid_step, 0.025)
    elseif scheme == "Schlafli"
        base_threshold = 0.3 * grid_step
        nnzz = CaseParams(base_threshold, 0.4 * grid_step, 0.05)
        nzzz = CaseParams(2.0 * base_threshold, 0.8 * grid_step, 0.025)
    else
        error("Unknown scheme: $scheme")
    end

    return WarpingSafetyParams(nnzz, nzzz)
end

"""
    slice_ambiguous_tetrahedra!(mesh::BlockMesh)

Slice tetrahedra crossing the isosurface (SDF zero level set).
Identifies elements crossing the boundary and replaces them with smaller
tetrahedra that accurately represent the surface.
"""
function slice_ambiguous_tetrahedra!(mesh::BlockMesh, scheme::String)
    @info "Slicing tetrahedra using trim_spikes logic..."

    warp_params = create_warping_params(scheme, mesh.grid_step)
    cut_map = Dict{Tuple{Int,Int},Int}()
    new_IEN = Vector{Vector{Int64}}()
    sizehint!(new_IEN, length(mesh.IEN))

    current_IEN = mesh.IEN
    mesh.IEN = Vector{Vector{Int64}}()

    for tet in current_IEN
        resulting_tets = apply_stencil_trim_spikes!(mesh, tet, cut_map, warp_params)
        for nt in resulting_tets
            push!(new_IEN, nt)
        end
    end

    mesh.IEN = new_IEN
    println("  After slicing: $(length(mesh.IEN)) tetrahedra")
end

# ----------------------------
# Helper function: Linear interpolation of the intersection with the zero level of SDF
# ----------------------------
function interpolate_zero(
    p1::SVector{3,Float64},
    p2::SVector{3,Float64},
    f1::Float64,
    f2::Float64,
    mesh::BlockMesh;
    tol = mesh.grid_tol,
    max_iter = 20,
)::Int64
    # coords of positive point, coords of negative point, sdf of positive point, sdf of negative point

    pos1 = f1 >= -tol
    pos2 = f2 >= -tol
    # If both points have the same polarity according to tolerance, interpolation cannot be done correctly.
    if pos1 == pos2
        error(
            "Both points have the same 'tolerance' polarity; one point must be close to zero (positive) and the other significantly negative.",
        )
    end

    # Initialize interval: low and high - assuming p1 and p2 are ordered by f
    low, high = p1, p2
    f_low, f_high = f1, f2
    mid = low
    for iter = 1:max_iter
        mid = (low + high) / 2.0
        f_mid = eval_sdf(mesh, mid)
        # If the value is close enough to zero, end the iteration
        if abs(f_mid) < tol
            break
        end
        # Update one of the interval endpoints according to the sign of f_mid
        if sign(f_mid) == sign(f_low)
            low, f_low = mid, f_mid
        else
            high, f_high = mid, f_mid
        end
    end

    # Quantize the found point to avoid duplicates in the hashtable
    p_key = quantize(mid, tol)
    sdf_of_iterp_point = eval_sdf(mesh, SVector{3,Float64}(p_key))
    # println("check interp sdf: ", sdf_of_iterp_point)
    if haskey(mesh.node_hash, p_key)
        return mesh.node_hash[p_key]
    else
        push!(mesh.X, mid)
        push!(mesh.node_sdf, 0.0)  # or set exactly to 0
        new_index = length(mesh.X)
        mesh.node_hash[p_key] = new_index
        return new_index
    end
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

Fix the orientation of an inverted tetrahedron by trying multiple vertex swap strategies.
Modifies the tetrahedron in-place. Returns true if fixed, false otherwise.

Tries three different strategies:
1. Swap vertices 3 and 4
2. Swap vertices 1 and 2
3. Swap vertices 2 and 3

Each strategy is tested and if successful, the fix is kept. If all strategies fail,
the tetrahedron is reverted to original state and a warning is issued.
"""
function fix_tetrahedron_orientation!(mesh::BlockMesh, tet::Vector{Int})
    if length(tet) != 4
        @warn "Attempting to fix orientation of non-tetrahedron: $tet"
        return false
    end

    # Check if already correctly oriented
    if check_tetrahedron_orientation(mesh, tet)
        return false  # No fix needed
    end

    # Store original ordering in case all strategies fail
    original_tet = copy(tet)

    # Strategy 1: Swap vertices 3 and 4
    tet[3], tet[4] = tet[4], tet[3]
    if check_tetrahedron_orientation(mesh, tet)
        return true  # Successfully fixed
    end

    # Revert strategy 1 and try strategy 2: Swap vertices 1 and 2
    tet[3], tet[4] = tet[4], tet[3]  # Revert to original
    tet[1], tet[2] = tet[2], tet[1]
    if check_tetrahedron_orientation(mesh, tet)
        println("strategy 1")
        return true  # Successfully fixed
    end

    # Revert strategy 2 and try strategy 3: Swap vertices 2 and 3
    tet[1], tet[2] = tet[2], tet[1]  # Revert to original
    tet[2], tet[3] = tet[3], tet[2]
    if check_tetrahedron_orientation(mesh, tet)
        println("strategy 2")
        return true  # Successfully fixed
    end

    # All strategies failed - revert to original and report failure
    tet[2], tet[3] = tet[3], tet[2]  # Revert to original
    @warn "Could not fix orientation for tetrahedron: $tet (vertices: $(mesh.X[tet[1]]), $(mesh.X[tet[2]]), $(mesh.X[tet[3]]), $(mesh.X[tet[4]]))"
    return false
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
        # @warn "Warp did not converge to surface: final SDF = $final_sdf"
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

# ----------------------------
# Case-specific warping handlers
# ----------------------------
"""
    process_nnzz_case!(mesh, s_idx, r_idx, q_idx, p_idx, sdf_s, sdf_r, cut_map, params, tol, flipped)

Process NNZZ case: two nodes outside (negative SDF), two on surface (zero SDF).

Warps both negative nodes to the surface with strict safety checks:
1. Original centroid must be near surface (threshold_distance)
2. Node displacements must be limited (max_node_displacement)
3. New volume must be reasonable (min_volume_ratio)
4. New centroid must be inside geometry (positive SDF)
5. Element must have correct orientation (positive Jacobian)

# Arguments
- `mesh::BlockMesh`: The tetrahedral mesh
- `s_idx, r_idx, q_idx, p_idx::Int`: Node indices (sorted by SDF: s ≤ r ≤ q ≤ p)
- `sdf_s, sdf_r::Float64`: SDF values at nodes s and r (negative)
- `cut_map::Dict`: Cache for avoiding duplicate node creation
- `params::CaseParams`: Safety parameters for this case
- `tol::Float64`: Tolerance for surface detection
- `flipped::Bool`: Indicates if element orientation was flipped during sorting

# Returns
- `Vector{Vector{Int64}}`: New tetrahedron if all checks pass, empty vector otherwise
"""
function process_nnzz_case!(
    mesh::BlockMesh,
    s_idx::Int,
    r_idx::Int,
    q_idx::Int,
    p_idx::Int,
    sdf_s::Float64,
    sdf_r::Float64,
    cut_map::Dict{Tuple{Int,Int},Int},
    params::CaseParams,
    tol::Float64,
    flipped::Bool,
)::Vector{Vector{Int64}}

    # Safety check 1: Original centroid distance from surface
    # Centrum původního elementu musí být blízko isosurface pro bezpečný warping
    original_centroid =
        (mesh.X[s_idx] + mesh.X[r_idx] + mesh.X[q_idx] + mesh.X[p_idx]) / 4.0
    if abs(eval_sdf(mesh, original_centroid)) > params.threshold_distance
        return Vector{Vector{Int64}}()
    end

    # Calculate original element volume for comparison
    # Potřebujeme pro kontrolu, že nový element nebude příliš malý
    v_orig = [mesh.X[s_idx], mesh.X[r_idx], mesh.X[q_idx], mesh.X[p_idx]]
    vol_orig =
        abs(
            dot(cross(v_orig[2] - v_orig[1], v_orig[3] - v_orig[1]), v_orig[4] - v_orig[1]),
        ) / 6.0

    # Warp both negative nodes to the isosurface
    # Oba uzly s negativním SDF se přesunou na isosurface
    s_warped = warp_node_to_surface(mesh, mesh.X[s_idx], sdf_s)
    r_warped = warp_node_to_surface(mesh, mesh.X[r_idx], sdf_r)

    # Safety check 2: Node displacement limits
    # Uzly se nesmí posunout příliš daleko (risk of mesh tangling)
    if norm(s_warped - mesh.X[s_idx]) > params.max_node_displacement ||
       norm(r_warped - mesh.X[r_idx]) > params.max_node_displacement
        return Vector{Vector{Int64}}()
    end

    # Add warped nodes to mesh (uses cache to avoid duplicates)
    # Přidá nové uzly do meshe, využívá hash pro eliminaci duplikátů
    s_new = add_warped_node!(mesh, s_warped, cut_map)
    r_new = add_warped_node!(mesh, r_warped, cut_map)

    # Create tetrahedron with correct orientation based on flipped flag
    # Respektuje původní orientaci elementu při vytváření nového tetraherdonu
    new_tet = if flipped
        [r_new, s_new, q_idx, p_idx]  # Flipped order - swap first two nodes
    else
        [s_new, r_new, q_idx, p_idx]  # Normal order
    end

    # Safety check 3: Volume ratio must be reasonable
    # Nový element nesmí být degenerovaný (příliš malý objem)
    v_new = [mesh.X[new_tet[1]], mesh.X[new_tet[2]], mesh.X[new_tet[3]], mesh.X[new_tet[4]]]
    vol_new =
        abs(dot(cross(v_new[2] - v_new[1], v_new[3] - v_new[1]), v_new[4] - v_new[1])) / 6.0

    if vol_new < params.min_volume_ratio * vol_orig
        return Vector{Vector{Int64}}()
    end

    # Safety check 4: New centroid must be inside geometry
    # Centrum nového elementu musí být uvnitř geometrie (pozitivní SDF)
    centroid =
        (
            mesh.X[new_tet[1]] +
            mesh.X[new_tet[2]] +
            mesh.X[new_tet[3]] +
            mesh.X[new_tet[4]]
        ) / 4.0
    if eval_sdf(mesh, centroid) <= tol
        return Vector{Vector{Int64}}()
    end

    # Safety check 5: Correct element orientation
    # Finální kontrola a případná oprava orientace (pozitivní Jacobián)
    fix_tetrahedron_orientation!(mesh, new_tet)

    # Verify orientation is correct after fix attempt
    if !check_tetrahedron_orientation(mesh, new_tet)
        @warn "NNZZ case: Failed to achieve positive orientation for tet $new_tet"
        return Vector{Vector{Int64}}()
    end

    return [new_tet]
end

"""
    process_nzzz_case!(mesh, s_idx, r_idx, q_idx, p_idx, sdf_s, cut_map, params, tol, flipped)

Process NZZZ case: one node outside (negative SDF), three on surface (zero SDF).

Warps single negative node to surface with relaxed safety checks (safer case).

# Arguments
- `mesh::BlockMesh`: The tetrahedral mesh
- `s_idx, r_idx, q_idx, p_idx::Int`: Node indices (sorted by SDF)
- `sdf_s::Float64`: SDF value at node s (negative)
- `cut_map::Dict`: Cache for avoiding duplicate node creation
- `params::CaseParams`: Safety parameters (relaxed for this case)
- `tol::Float64`: Tolerance for surface detection
- `flipped::Bool`: Indicates if element orientation was flipped during sorting

# Returns
- `Vector{Vector{Int64}}`: New tetrahedron if all checks pass, empty vector otherwise
"""
function process_nzzz_case!(
    mesh::BlockMesh,
    s_idx::Int,
    r_idx::Int,
    q_idx::Int,
    p_idx::Int,
    sdf_s::Float64,
    cut_map::Dict{Tuple{Int,Int},Int},
    params::CaseParams,
    tol::Float64,
    flipped::Bool,
)::Vector{Vector{Int64}}

    # Safety check 1: Original centroid distance from surface (relaxed threshold)
    # Relaxed threshold protože warpujeme jen jeden uzel (bezpečnější)
    original_centroid =
        (mesh.X[s_idx] + mesh.X[r_idx] + mesh.X[q_idx] + mesh.X[p_idx]) / 4.0
    if abs(eval_sdf(mesh, original_centroid)) > params.threshold_distance
        return Vector{Vector{Int64}}()
    end

    # Calculate original element volume for comparison
    v_orig = [mesh.X[s_idx], mesh.X[r_idx], mesh.X[q_idx], mesh.X[p_idx]]
    vol_orig =
        abs(
            dot(cross(v_orig[2] - v_orig[1], v_orig[3] - v_orig[1]), v_orig[4] - v_orig[1]),
        ) / 6.0

    # Warp single negative node to the isosurface
    # Pouze jeden uzel s negativním SDF se přesune na isosurface
    s_warped = warp_node_to_surface(mesh, mesh.X[s_idx], sdf_s)

    # Safety check 2: Node displacement limit (relaxed for single node)
    # Relaxed limit protože warpujeme jen jeden uzel
    if norm(s_warped - mesh.X[s_idx]) > params.max_node_displacement
        return Vector{Vector{Int64}}()
    end

    # Add warped node to mesh
    s_new = add_warped_node!(mesh, s_warped, cut_map)

    # Create tetrahedron with correct orientation based on flipped flag
    # Pro NZZZ case: swap s_new s r_idx když je flipped
    new_tet = if flipped
        [r_idx, s_new, q_idx, p_idx]  # Flipped order - swap first two nodes
    else
        [s_new, r_idx, q_idx, p_idx]  # Normal order
    end

    # println(flipped)

    # Safety check 3: Volume ratio (relaxed minimum)
    # Relaxed threshold protože deformace je menší (jen jeden warpovaný uzel)
    v_new = [mesh.X[new_tet[1]], mesh.X[new_tet[2]], mesh.X[new_tet[3]], mesh.X[new_tet[4]]]
    vol_new =
        abs(dot(cross(v_new[2] - v_new[1], v_new[3] - v_new[1]), v_new[4] - v_new[1])) / 6.0

    if vol_new < params.min_volume_ratio * vol_orig
        return Vector{Vector{Int64}}()
    end

    # Safety check 4: New centroid must be inside geometry
    # Centrum musí zůstat uvnitř geometrie
    centroid =
        (
            mesh.X[new_tet[1]] +
            mesh.X[new_tet[2]] +
            mesh.X[new_tet[3]] +
            mesh.X[new_tet[4]]
        ) / 4.0
    if eval_sdf(mesh, centroid) <= tol
        return Vector{Vector{Int64}}()
    end

    # Safety check 5: Correct element orientation
    # Finální kontrola orientace s možností opravy
    fix_tetrahedron_orientation!(mesh, new_tet)

    # Verify orientation is correct after fix attempt
    if !check_tetrahedron_orientation(mesh, new_tet)
        @warn "NZZZ case: Failed to achieve positive orientation for tet $new_tet"
        return Vector{Vector{Int64}}()
    end

    return [new_tet]
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
    warp_params::WarpingSafetyParams,
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

    # --- Case analysis based on SDF sign patterns --- (S ≤ R ≤ Q ≤ P)
    if is_s_neg && is_p_zero # Cases: NNNZ, NNZZ, NZZZ
        if is_s_neg && is_r_neg && is_q_neg && is_p_zero
            # NNNZ: Three nodes outside, one on surface
            return Vector{Vector{Int64}}()

        elseif is_s_neg && is_r_neg && is_q_zero && is_p_zero
            # NNZZ: Two nodes outside, two on surface
            return process_nnzz_case!(
                mesh,
                s_idx,
                r_idx,
                q_idx,
                p_idx,
                sdf_s,
                sdf_r,
                cut_map,
                warp_params.nnzz,
                tol,
                flipped,  # ← PŘIDÁN PARAMETR
            )

        elseif is_s_neg && is_r_zero && is_q_zero && is_p_zero
            # NZZZ: One node outside, three on surface
            return process_nzzz_case!(
                mesh,
                s_idx,
                r_idx,
                q_idx,
                p_idx,
                sdf_s,
                cut_map,
                warp_params.nzzz,
                tol,
                flipped,  # ← PŘIDÁN PARAMETR
            )
        else
            println("  -> NO MATCH!")
            return Vector{Vector{Int64}}()
        end

        # Case (Z,!N,!N,P): Surface node with all others on surface or inside
    elseif is_s_zero && !is_r_neg && !is_q_neg && is_p_pos # ZZZP, ZZPP, ZPPP
        # println("ZZZP, ZZPP, ZPPP")
        return [tet]

        # Case NNNP: Three nodes outside, one inside
    elseif is_s_neg && is_r_neg && is_q_neg && is_p_pos
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
