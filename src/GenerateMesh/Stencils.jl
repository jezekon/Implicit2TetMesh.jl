# Functions for tetrahedral mesh slicing along an isosurface
# Implements the "trim spikes" algorithm for accurate surface representation

"""
    remove_inverted_elements!(mesh::BlockMesh)

Removes tetrahedral elements with negative Jacobian determinants (negative volumes)
from the mesh. These inverted elements can occur during the warping process 
when nodes are adjusted to fit the isosurface.

# Arguments
- `mesh::BlockMesh`: The mesh to be cleaned

# Returns
- `mesh::BlockMesh`: The modified mesh with inverted elements removed

# Note
This function should be called before updating the mesh connectivity to ensure
that all inverted elements are properly removed from the mesh.
"""
function remove_inverted_elements!(mesh::BlockMesh)
    # Set default tolerance if not provided
    # if volume_tolerance === nothing
        volume_tolerance = mesh.grid_tol * 1e-6
    # end
    
    # Pre-allocate output for better performance
    valid_elements = Vector{Vector{Int}}()
    sizehint!(valid_elements, length(mesh.IEN))
    
    # Track problematic elements for reporting
    negative_elements = Tuple{Int,Float64}[]
    zero_elements = Tuple{Int,Float64}[]
    
    for (elem_idx, tet) in enumerate(mesh.IEN)
        # Get vertices of the tetrahedron as static vectors for better performance
        vertices = SVector{4,SVector{3,Float64}}(
            mesh.X[tet[1]], mesh.X[tet[2]], mesh.X[tet[3]], mesh.X[tet[4]]
        )
        
        # Calculate volume using cross product method
        a = vertices[2] - vertices[1]
        b = vertices[3] - vertices[1]
        c = vertices[4] - vertices[1]
        
        # Calculate determinant (6 times the volume)
        det_value = dot(a, cross(b, c))
        
        # Actual volume
        volume = det_value / 6.0
        
        # Categorize element based on volume
        if det_value < 0
            # Negative volume (inverted element)
            push!(negative_elements, (elem_idx, volume))
        elseif abs(det_value) <= volume_tolerance
            # Near-zero volume (degenerate element)
            push!(zero_elements, (elem_idx, volume))
        else
            # Valid element with positive volume
            push!(valid_elements, tet)
        end
    end
    
    # Sort problematic elements by volume for reporting
    sort!(negative_elements, by=x->x[2])
    sort!(zero_elements, by=x->abs(x[2]))
    
    # Get total count of removed elements
    negative_count = length(negative_elements)
    zero_count = length(zero_elements)
    total_removed = negative_count + zero_count
    
    # Report problematic elements if any were found
    if total_removed > 0
        @info "Removing $total_removed problematic elements from the mesh ($negative_count negative, $zero_count zero-volume)."
        
        # Report negative elements
        if !isempty(negative_elements)
            println("⚠️  REMOVING NEGATIVE VOLUME ELEMENTS ⚠️")
            println("+------------+------------------------+")
            println("| Element ID | Volume                 |")
            println("+------------+------------------------+")
            
            # Show up to 5 most negative elements
            for (i, (elem_id, volume)) in enumerate(negative_elements[1:min(5, length(negative_elements))])
                println(@sprintf("| %-10d | %-22.6e |", elem_id, volume))
            end
            
            if length(negative_elements) > 5
                println("| ... and $(length(negative_elements) - 5) more negative volume elements")
            end
            println("+------------+------------------------+")
        end
        
        # Report zero volume elements
        if !isempty(zero_elements)
            println("⚠️  REMOVING ZERO VOLUME ELEMENTS ⚠️")
            println("+------------+------------------------+")
            println("| Element ID | Volume                 |")
            println("+------------+------------------------+")
            
            # Show up to 5 elements with smallest absolute volume
            for (i, (elem_id, volume)) in enumerate(zero_elements[1:min(5, length(zero_elements))])
                println(@sprintf("| %-10d | %-22.6e |", elem_id, volume))
            end
            
            if length(zero_elements) > 5
                println("| ... and $(length(zero_elements) - 5) more zero volume elements")
            end
            println("+------------+------------------------+")
        end
    else
        @info "No problematic elements found in the mesh."
    end
    
    # Update mesh connectivity
    mesh.IEN = valid_elements
    
    return mesh
end

"""
    slice_ambiguous_tetrahedra!(mesh::BlockMesh)

Slice tetrahedra crossing the isosurface (SDF zero level set).
Identifies elements crossing the boundary and replaces them with smaller
tetrahedra that accurately represent the surface.
"""
function slice_ambiguous_tetrahedra!(mesh::BlockMesh)
    @info "Slicing tetrahedra using trim_spikes logic..."
    new_IEN = Vector{Vector{Int64}}()
    sizehint!(new_IEN, length(mesh.IEN)) # Pre-allocate for performance

    # Cache for storing cut edge results to avoid duplicate cutting
    cut_map = Dict{Tuple{Int, Int}, Int}()

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
    new_node_count = length(mesh.X)
    new_tet_count = length(mesh.IEN)

    @info "After trim_spikes slicing: $(new_tet_count) tetrahedra (added $(new_tet_count - original_tet_count)), $(new_node_count) nodes (added $(new_node_count - original_node_count))"

    # Note: After slicing, update_connectivity!(mesh) should be called externally
end

"""
    cut_edge!(i::Int, j::Int, mesh::BlockMesh, node_sdf::Vector{Float64}, cut_map::Dict{Tuple{Int, Int}, Int})::Int

Calculate the intersection point between an edge (i,j) and the isosurface.
Returns the node index at the intersection point, creating a new node if needed.
Handles special cases where nodes are already on the surface.
"""
function cut_edge!(
    i::Int, j::Int,
    mesh::BlockMesh,
    node_sdf::Vector{Float64},
    cut_map::Dict{Tuple{Int, Int}, Int}
)::Int
    # Get SDF values for both endpoints
    sdf_i = node_sdf[i]
    sdf_j = node_sdf[j]

    tol = mesh.grid_tol
    is_i_zero = abs(sdf_i) < tol
    is_j_zero = abs(sdf_j) < tol

    # Handle cases where nodes are already on the surface
    if is_i_zero && is_j_zero
        # Both on surface - return the lower index for consistency
        return min(i, j)
    elseif is_i_zero
        return i # Node i is already on the surface
    elseif is_j_zero
        return j # Node j is already on the surface
    end

    # Ensure nodes are on opposite sides of the isosurface
    if sign(sdf_i) == sign(sdf_j)
        error("cut_edge! called with nodes on same side of isosurface: i=$i (sdf=$(sdf_i)), j=$j (sdf=$(sdf_j))")
    end

    # Canonical edge representation for consistent map lookup
    edge = (min(i, j), max(i, j))

    # Check if this edge has already been cut
    if haskey(cut_map, edge)
        return cut_map[edge]
    end

    # Get node positions
    pos_i = mesh.X[i]
    pos_j = mesh.X[j]

    # Linear interpolation to find intersection point
    alpha = sdf_i / (sdf_i - sdf_j)
    alpha = clamp(alpha, 0.0, 1.0) # Handle numerical imprecision

    # Calculate intersection position
    new_pos = (1.0 - alpha) * pos_i + alpha * pos_j

    # Quantize position for consistent node identification
    p_key = quantize(new_pos, mesh.grid_tol)

    # Either use existing node or create new one
    local new_index::Int
    if haskey(mesh.node_hash, p_key)
        # Reuse existing node
        new_index = mesh.node_hash[p_key]
        # Ensure its SDF is exactly zero if it's a surface node
        if abs(mesh.node_sdf[new_index]) > tol
             mesh.node_sdf[new_index] = 0.0
        end
    else
        # Create new node at the intersection
        push!(mesh.X, new_pos)
        push!(mesh.node_sdf, 0.0) # Surface node (SDF = 0)
        new_index = length(mesh.X)
        mesh.node_hash[p_key] = new_index
    end

    # Cache result for future edge cuts
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
    vertices = [mesh.X[tet[i]] for i in 1:4]

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
    apply_stencil_trim_spikes!(mesh::BlockMesh, tet::Vector{Int64}, cut_map::Dict{Tuple{Int, Int}, Int})

Apply trimming algorithm to a tetrahedron that crosses the isosurface.
Returns new tetrahedra that accurately represent the interior region (SDF ≥ 0).
Uses case-by-case analysis based on SDF sign patterns for robust slicing.
"""
function apply_stencil_trim_spikes!(
    mesh::BlockMesh,
    tet::Vector{Int64},
    cut_map::Dict{Tuple{Int, Int}, Int}
)::Vector{Vector{Int64}}

    # Get SDF values for tetrahedron nodes
    node_indices = tet
    node_sdf = [mesh.node_sdf[idx] for idx in node_indices]
    tol = mesh.grid_tol

    # --- Special cases handling ---
    # Case 1: All nodes outside (SDF < -tol) - discard tetrahedron
    if all(s -> s < -tol, node_sdf)
        return Vector{Vector{Int64}}()
    end
    
    # Case 2: All nodes inside or on surface (SDF ≥ -tol) - keep original
    if all(s -> s >= -tol, node_sdf)
        return [tet]
    end
    
    # --- General case: Tetrahedron crossing the isosurface ---
    # Sort vertices by SDF value (smallest first)
    vert_data = [(node_sdf[i], node_indices[i]) for i in 1:4]
    p = [1, 2, 3, 4]
    flipped = false
    
    # Comparison function for consistent vertex ordering
    less_than(idx1, idx2) = (vert_data[idx1][1] < vert_data[idx2][1]) || 
                            (vert_data[idx1][1] == vert_data[idx2][1] && 
                             vert_data[idx1][2] < vert_data[idx2][2])

    # Sort vertices while tracking orientation flips
    if less_than(p[2], p[1]) p[1], p[2] = p[2], p[1]; flipped = !flipped end
    if less_than(p[4], p[3]) p[3], p[4] = p[4], p[3]; flipped = !flipped end
    if less_than(p[3], p[1]) p[1], p[3] = p[3], p[1]; flipped = !flipped end
    if less_than(p[4], p[2]) p[2], p[4] = p[4], p[2]; flipped = !flipped end
    if less_than(p[3], p[2]) p[2], p[3] = p[3], p[2]; flipped = !flipped end

    # Sorted vertices (s has lowest SDF, p has highest)
    s_idx, r_idx, q_idx, p_idx = [vert_data[pi][2] for pi in p]
    sdf_s, sdf_r, sdf_q, sdf_p = [vert_data[pi][1] for pi in p]

    # Classify vertices by SDF sign (N=negative, P=positive, Z=zero)
    is_s_neg = sdf_s < -tol
    is_r_neg = sdf_r < -tol
    is_q_neg = sdf_q < -tol

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

    # --- Case analysis based on SDF sign patterns ---
    
    # Case ZPPPP: Surface node with all others inside - special handling
    if is_s_zero && !is_r_neg && !is_q_neg && !is_p_neg
        # Subcase ZZZZ: All nodes on surface, check centroid
        if is_r_zero && is_q_zero && is_p_zero
             centroid = (mesh.X[s_idx] + mesh.X[r_idx] + mesh.X[q_idx] + mesh.X[p_idx]) / 4.0
             if eval_sdf(mesh, centroid) < -tol
                  return Vector{Vector{Int64}}() # Discard if centroid outside
             else
                  return Vector{Vector{Int64}}() # Discard surface tetrahedra
             end
        else # Other Z+++ patterns - discard
             return Vector{Vector{Int64}}()
        end

    # Case NNNP: Three nodes outside, one inside
    elseif is_s_neg && is_r_neg && is_q_neg && (is_p_pos || is_p_zero)
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
    elseif is_s_neg && (is_p_pos || is_p_zero) && (is_q_pos || is_q_zero) && (is_r_pos || is_r_zero)
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
    elseif is_s_neg && is_r_neg && (is_q_pos || is_q_zero) && (is_p_pos || is_p_zero)
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
    elseif is_s_neg && is_r_zero && (is_q_pos || is_q_zero) && (is_p_pos || is_p_zero)
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
    elseif is_s_neg && is_r_neg && is_q_zero && (is_p_pos || is_p_zero)
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
    elseif is_s_neg && is_r_zero && is_q_zero && (is_p_pos || is_p_zero)
        # Cut edge from outside node to inside node
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        
        # Create one tetrahedron
        if flipped
             add_tet!([q_idx, r_idx, p_idx, sp])
        else
             add_tet!([r_idx, q_idx, p_idx, sp])
        end

    # Unexpected SDF pattern
    else
        @warn "Unexpected SDF pattern in apply_stencil_trim_spikes!: s=$sdf_s, r=$sdf_r, q=$sdf_q, p=$sdf_p. Indices: s=$s_idx, r=$r_idx, q=$q_idx, p=$p_idx. Tet: $tet. Flipped: $flipped"
        return Vector{Vector{Int64}}() # Discard in case of unexpected pattern
    end

    return new_tets
end
