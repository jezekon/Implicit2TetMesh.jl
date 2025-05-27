"""
    calculate_mesh_volume(mesh::BlockMesh) -> Float64

Vypočítá celkový objem tetrahedrální sítě jako sumu objemů jednotlivých tetrahedrů.
"""
function calculate_mesh_volume(mesh::BlockMesh)
    total_volume = 0.0
    
    for tet in mesh.IEN
        # Získej vrcholy tetrahedru
        v1 = mesh.X[tet[1]]
        v2 = mesh.X[tet[2]]
        v3 = mesh.X[tet[3]]
        v4 = mesh.X[tet[4]]
        
        # Vypočítaj objem tetrahedru: |det(v2-v1, v3-v1, v4-v1)| / 6
        edge1 = v2 - v1
        edge2 = v3 - v1
        edge3 = v4 - v1
        
        volume = abs(dot(edge1, cross(edge2, edge3))) / 6.0
        total_volume += volume
    end
    
    return total_volume
end

"""
    detect_surface_nodes(mesh::BlockMesh, sdf_tolerance::Float64=0.01) -> Set{Int}

Detects nodes that are actually on or near the isosurface (SDF ≈ 0),
rather than just boundary nodes of the mesh topology.
This is more appropriate for volume correction.
"""
function detect_surface_nodes(mesh::BlockMesh, sdf_tolerance::Float64=0.01)
    surface_nodes = Set{Int}()
    
    # Find nodes with SDF values close to zero (on the isosurface)
    for (node_idx, sdf_val) in enumerate(mesh.node_sdf)
        if abs(sdf_val) <= sdf_tolerance
            push!(surface_nodes, node_idx)
        end
    end
    
    @info "Detected $(length(surface_nodes)) surface nodes (SDF ≈ 0)"
    return surface_nodes
end

"""
    check_element_validity(mesh::BlockMesh) -> (Int, Int)

Checks for inverted or degenerate elements in the mesh.
Returns (negative_count, zero_volume_count).
"""
function check_element_validity(mesh::BlockMesh)
    negative_count = 0
    zero_volume_count = 0
    tolerance = 1e-12
    
    for tet in mesh.IEN
        # Skip degenerate elements with duplicate vertices
        if length(Set(tet)) != 4
            zero_volume_count += 1
            continue
        end
        
        # Calculate determinant (proportional to signed volume)
        vertices = [mesh.X[node_idx] for node_idx in tet]
        a = vertices[2] - vertices[1]
        b = vertices[3] - vertices[1]
        c = vertices[4] - vertices[1]
        det_value = dot(a, cross(b, c))
        
        if det_value < -tolerance
            negative_count += 1
        elseif abs(det_value) <= tolerance
            zero_volume_count += 1
        end
    end
    
    return negative_count, zero_volume_count
end

"""
    correct_mesh_volume!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3}; 
                               tolerance::Float64=0.005, max_iterations::Int=20)

Improved volume correction method that:
1. Uses surface nodes (SDF ≈ 0) instead of boundary nodes
2. Applies smaller, more conservative adjustments
3. Validates mesh integrity after each adjustment
4. Uses a more robust bisection approach
"""
function correct_mesh_volume!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3}; 
                                   tolerance::Float64=0.005, max_iterations::Int=20)
    @info "Starting improved volume correction..."
    
    # Calculate reference volume from SDF data
    reference_volume = Float64(calculate_volume_from_sdf(fine_sdf, fine_grid))
    @info "Reference volume from SDF: $reference_volume"
    
    # Detect surface nodes (nodes with SDF ≈ 0)
    surface_nodes = detect_surface_nodes(mesh, 0.02)  # Slightly larger tolerance
    
    if isempty(surface_nodes)
        @warn "No surface nodes detected. Cannot perform volume correction."
        return false
    end
    
    # Store original positions of surface nodes
    original_positions = Dict{Int, SVector{3,Float64}}()
    original_sdf_values = Dict{Int, Float64}()
    for node_idx in surface_nodes
        original_positions[node_idx] = mesh.X[node_idx]
        original_sdf_values[node_idx] = mesh.node_sdf[node_idx]
    end
    
    # Calculate initial mesh volume and error
    initial_volume = calculate_mesh_volume(mesh)
    initial_error = abs(initial_volume - reference_volume) / reference_volume
    @info "Initial mesh volume: $initial_volume"
    @info "Initial relative error: $(initial_error * 100)%"
    
    if initial_error <= tolerance
        @info "Mesh volume already within tolerance"
        return true
    end
    
    # Check initial mesh validity
    neg_count, zero_count = check_element_validity(mesh)
    @info "Initial mesh validity: $neg_count negative elements, $zero_count zero-volume elements"
    
    # Determine adjustment direction and magnitude
    volume_deficit = reference_volume - initial_volume
    adjustment_direction = sign(volume_deficit)  # +1 for expanding, -1 for contracting
    
    @info "Volume deficit: $volume_deficit ($(adjustment_direction > 0 ? "need to expand" : "need to contract"))"
    
    # Use smaller, more conservative steps
    max_sdf_adjustment = 0.02  # Much smaller than before
    current_adjustment = 0.0
    step_size = max_sdf_adjustment / max_iterations
    
    best_volume = initial_volume
    best_error = initial_error
    best_adjustment = 0.0
    
    # Progressive adjustment approach
    for iteration in 1:max_iterations
        current_adjustment = step_size * iteration * adjustment_direction
        
        @info "Iteration $iteration: Testing SDF adjustment $current_adjustment"
        
        # Restore original positions
        restore_surface_nodes!(mesh, surface_nodes, original_positions, original_sdf_values)
        
        # Apply conservative adjustment to surface nodes
        successful_moves = conservative_surface_adjustment!(mesh, surface_nodes, current_adjustment)
        
        if successful_moves == 0
            @warn "No nodes could be adjusted at level $current_adjustment"
            break
        end
        
        # Check mesh validity after adjustment
        neg_count, zero_count = check_element_validity(mesh)
        if neg_count > 0
            @warn "Adjustment created $neg_count inverted elements, reverting"
            restore_surface_nodes!(mesh, surface_nodes, original_positions, original_sdf_values)
            break
        end
        
        # Calculate new volume and error
        new_volume = calculate_mesh_volume(mesh)
        new_error = abs(new_volume - reference_volume) / reference_volume
        
        @info "  New volume: $new_volume"
        @info "  Relative error: $(new_error * 100)%"
        @info "  Successfully adjusted: $successful_moves nodes"
        
        # Track best solution
        if new_error < best_error
            best_error = new_error
            best_volume = new_volume
            best_adjustment = current_adjustment
        end
        
        # Check for convergence
        if new_error <= tolerance
            @info "Volume correction converged after $iteration iterations"
            return true
        end
        
        # If error is getting worse, try smaller steps
        if new_error > initial_error * 1.1  # 10% worse than initial
            @warn "Error increasing significantly, reducing step size"
            step_size *= 0.5
            if step_size < 1e-6
                break
            end
        end
    end
    
    # Apply best solution found
    @info "Applying best solution with adjustment $best_adjustment (error: $(best_error * 100)%)"
    restore_surface_nodes!(mesh, surface_nodes, original_positions, original_sdf_values)
    conservative_surface_adjustment!(mesh, surface_nodes, best_adjustment)
    
    # Final assessment
    final_volume = calculate_mesh_volume(mesh)
    final_error = abs(final_volume - reference_volume) / reference_volume
    @info "Final volume: $final_volume"
    @info "Final relative error: $(final_error * 100)%"
    
    return final_error <= tolerance
end

"""
    restore_surface_nodes!(mesh::BlockMesh, surface_nodes::Set{Int}, 
                          original_positions::Dict{Int, SVector{3,Float64}},
                          original_sdf_values::Dict{Int, Float64})

Restores surface nodes to their original positions and SDF values.
"""
function restore_surface_nodes!(mesh::BlockMesh, surface_nodes::Set{Int}, 
                                original_positions::Dict{Int, SVector{3,Float64}},
                                original_sdf_values::Dict{Int, Float64})
    for node_idx in surface_nodes
        mesh.X[node_idx] = original_positions[node_idx]
        mesh.node_sdf[node_idx] = original_sdf_values[node_idx]
    end
end

"""
    conservative_surface_adjustment!(mesh::BlockMesh, surface_nodes::Set{Int}, 
                                    sdf_adjustment::Float64; max_displacement::Float64=0.01) -> Int

Applies conservative adjustment to surface nodes.
Returns the number of successfully adjusted nodes.
"""
function conservative_surface_adjustment!(mesh::BlockMesh, surface_nodes::Set{Int}, 
                                        sdf_adjustment::Float64; max_displacement::Float64=0.01)
    nodes_adjusted = 0
    
    for node_idx in surface_nodes
        # Calculate gradient direction for this node
        gradient = compute_gradient(mesh, mesh.X[node_idx])
        gradient_norm = norm(gradient)
        
        if gradient_norm < 1e-10
            continue  # Skip nodes with zero gradient
        end
        
        # Normalize gradient
        normalized_gradient = gradient / gradient_norm
        
        # Calculate displacement based on desired SDF change
        # For positive sdf_adjustment (expanding), move outward (in gradient direction)
        # For negative sdf_adjustment (contracting), move inward (opposite to gradient)
        displacement_magnitude = abs(sdf_adjustment) / gradient_norm
        displacement_magnitude = min(displacement_magnitude, max_displacement)
        
        displacement_vector = sign(sdf_adjustment) * displacement_magnitude * normalized_gradient
        
        # Apply displacement
        new_position = mesh.X[node_idx] + displacement_vector
        
        # Update node position and recalculate SDF
        mesh.X[node_idx] = new_position
        mesh.node_sdf[node_idx] = eval_sdf(mesh, new_position)
        
        nodes_adjusted += 1
    end
    
    return nodes_adjusted
end

"""
    alternative_volume_correction!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3}; 
                                  tolerance::Float64=0.005)

Alternative approach: uniform scaling of the entire mesh.
This is more predictable but less sophisticated.
"""
function alternative_volume_correction!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3}; 
                                      tolerance::Float64=0.005)
    @info "Starting alternative volume correction (uniform scaling)..."
    
    reference_volume = Float64(calculate_volume_from_sdf(fine_sdf, fine_grid))
    current_volume = calculate_mesh_volume(mesh)
    
    @info "Reference volume: $reference_volume"
    @info "Current volume: $current_volume"
    
    # Calculate required scaling factor
    volume_ratio = reference_volume / current_volume
    scale_factor = cbrt(volume_ratio)  # Cube root for 3D scaling
    
    @info "Required scale factor: $scale_factor"
    
    # Find mesh centroid
    centroid = sum(mesh.X) / length(mesh.X)
    @info "Mesh centroid: $centroid"
    
    # Apply uniform scaling about centroid
    for i in 1:length(mesh.X)
        # Vector from centroid to node
        offset = mesh.X[i] - centroid
        # Scale the offset
        scaled_offset = scale_factor * offset
        # Update node position
        mesh.X[i] = centroid + scaled_offset
        # Recalculate SDF
        mesh.node_sdf[i] = eval_sdf(mesh, mesh.X[i])
    end
    
    # Check final volume
    final_volume = calculate_mesh_volume(mesh)
    final_error = abs(final_volume - reference_volume) / reference_volume
    
    @info "Final volume after scaling: $final_volume"
    @info "Final relative error: $(final_error * 100)%"
    
    return final_error <= tolerance
end

"""
    assess_volume_accuracy(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3})

Vyhodnotí přesnost objemu sítě ve srovnání s referenčním objemem ze SDF dat.
"""
function assess_volume_accuracy(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3})
    reference_volume = Float64(calculate_volume_from_sdf(fine_sdf, fine_grid))
    mesh_volume = calculate_mesh_volume(mesh)
    
    absolute_error = abs(mesh_volume - reference_volume)
    relative_error = absolute_error / reference_volume
    
    println("Volume Assessment:")
    println("  Reference volume (SDF): $reference_volume")
    println("  Mesh volume: $mesh_volume")
    println("  Absolute error: $absolute_error")
    println("  Relative error: $(relative_error * 100)%")
    
    return (reference_volume, mesh_volume, relative_error)
end
