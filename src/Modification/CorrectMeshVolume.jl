"""
    correct_mesh_volume!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3};
                        tolerance::Float64=0.01, max_iterations::Int=20)

Corrects the tetrahedral mesh volume to match the implicit geometry volume by adjusting 
surface nodes along SDF gradient directions using bisection method.

# Arguments
- `mesh::BlockMesh`: The tetrahedral mesh to be corrected
- `fine_sdf::Array{Float32,3}`: SDF values on fine grid for reference volume calculation
- `fine_grid::Array{Vector{Float32},3}`: Fine grid coordinates for reference volume calculation
- `tolerance::Float64=0.01`: Relative volume error tolerance for convergence (1% default)
- `max_iterations::Int=20`: Maximum number of bisection iterations

# Algorithm
1. Calculate reference volume from SDF data using calculate_volume_from_sdf
2. Detect surface nodes (nodes with SDF ≈ 0)
3. Determine initial adjustment direction: ±0.05 SDF level based on volume deficit
4. Use bisection method to find optimal SDF level that matches reference volume
5. Move surface nodes along SDF gradient (normal) direction to target SDF level

# Output
Prints volume and SDF level information for each iteration
"""
function correct_mesh_volume!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3};
                              tolerance::Float64=0.0001, max_iterations::Int=20)
    
    # Step 1: Calculate reference volume from SDF data
    reference_volume = Float64(calculate_volume_from_sdf(fine_sdf, fine_grid))
    @info "Reference volume from SDF data: $reference_volume"
    
    # Step 2: Detect surface nodes (nodes with SDF ≈ 0)
    surface_nodes = Int[]
    surface_tolerance = 0.02  # Tolerance for surface node detection
    
    for i in 1:length(mesh.X)
        if abs(mesh.node_sdf[i]) <= surface_tolerance
            push!(surface_nodes, i)
        end
    end
    
    @info "Detected $(length(surface_nodes)) surface nodes"
    
    if isempty(surface_nodes)
        @warn "No surface nodes found - cannot perform volume correction"
        return false
    end
    
    # Step 3: Helper function to calculate current mesh volume
    function calculate_current_mesh_volume()
        total_volume = 0.0
        for tet in mesh.IEN
            # Get tetrahedron vertices
            vertices = [mesh.X[tet[i]] for i in 1:4]
            
            # Calculate tetrahedron volume: |det(v2-v1, v3-v1, v4-v1)| / 6
            edge1 = vertices[2] - vertices[1]
            edge2 = vertices[3] - vertices[1] 
            edge3 = vertices[4] - vertices[1]
            volume = abs(dot(edge1, cross(edge2, edge3))) / 6.0
            total_volume += volume
        end
        return total_volume
    end
    
    # Step 4: Helper function to move surface nodes to target SDF level
    function move_surface_nodes_to_sdf_level!(target_sdf::Float64)
        for node_idx in surface_nodes
            # Calculate SDF gradient (surface normal) at current node position
            gradient = compute_gradient(mesh, mesh.X[node_idx])
            grad_norm = norm(gradient)
            
            # Skip nodes with zero gradient (degenerate cases)
            if grad_norm < 1e-10
                continue
            end
            
            # Calculate current SDF value at node
            current_sdf = eval_sdf(mesh, mesh.X[node_idx])
            
            # Calculate required SDF change to reach target level
            sdf_change = target_sdf - current_sdf
            
            # Calculate displacement along gradient direction
            # Positive SDF change = move outward (along gradient)
            # Negative SDF change = move inward (against gradient)
            unit_normal = gradient / grad_norm
            displacement = (sdf_change / grad_norm) * unit_normal
            
            # Update node position
            mesh.X[node_idx] += displacement
            
            # Update node SDF value
            mesh.node_sdf[node_idx] = eval_sdf(mesh, mesh.X[node_idx])
        end
    end
    
    # Step 5: Determine initial volume deficit and correction direction
    current_volume = calculate_current_mesh_volume()
    volume_deficit = reference_volume - current_volume
    
    @info "Initial mesh volume: $current_volume"
    @info "Volume deficit: $volume_deficit"
    
    # Determine initial SDF target based on volume deficit
    initial_sdf = volume_deficit > 0 ? -0.05 : 0.05
    correction_direction = volume_deficit > 0 ? "expand inward" : "contract outward"
    
    @info "Volume correction needed: $correction_direction"
    @info "Initial SDF target: $initial_sdf"
    
    # Step 6: Set up bisection interval
    if volume_deficit > 0
        # Need to expand volume: move inward (negative SDF direction)
        low_sdf, high_sdf = initial_sdf, 0.0
    else
        # Need to contract volume: move outward (positive SDF direction)  
        low_sdf, high_sdf = 0.0, initial_sdf
    end
    
    # Backup original node positions and SDF values for restoration
    original_positions = [mesh.X[i] for i in surface_nodes]
    original_sdf_values = [mesh.node_sdf[i] for i in surface_nodes]
    
    # Track best solution found
    best_volume = current_volume
    best_sdf = 0.0
    best_error = abs(current_volume - reference_volume) / reference_volume
    
    # Step 7: Bisection iterations to find optimal SDF level
    println("\nStarting bisection iterations:")
    println("Iter | SDF Level | Mesh Volume | Rel. Error")
    println("-----|-----------|-------------|------------")
    
    for iter in 1:max_iterations
        # Calculate midpoint SDF level
        mid_sdf = (low_sdf + high_sdf) / 2.0
        
        # Restore original node positions before each test
        for (idx, node_idx) in enumerate(surface_nodes)
            mesh.X[node_idx] = original_positions[idx]
            mesh.node_sdf[node_idx] = original_sdf_values[idx]
        end
        
        # Move nodes to current test SDF level
        move_surface_nodes_to_sdf_level!(mid_sdf)
        
        # Calculate resulting mesh volume
        new_volume = calculate_current_mesh_volume()
        relative_error = abs(new_volume - reference_volume) / reference_volume
        
        # Print iteration results 
        println("$(lpad(iter,4)) | $(round(mid_sdf, digits=6)) | $(round(new_volume, digits=6)) | $(round(relative_error * 100, digits=2))%")
        
        # Track best solution encountered
        if relative_error < best_error
            best_volume = new_volume
            best_sdf = mid_sdf
            best_error = relative_error
        end
        
        # Check for convergence
        if relative_error <= tolerance
            @info "Volume correction converged after $iter iterations!"
            return true
        end
        
        # Update bisection interval based on volume comparison
        if new_volume > reference_volume
           low_sdf = mid_sdf  # Move toward more negative SDF (more inward)
        else
           high_sdf = mid_sdf # Move toward less negative SDF (less inward)
        end
        
        # Check if interval became too small (convergence unlikely)
        if abs(high_sdf - low_sdf) < 1e-8
            @warn "SDF interval too small ($((high_sdf - low_sdf))) - stopping iterations"
            break
        end
    end
    
    # Step 8: Apply best solution found
    @info "Applying best solution found:"
    @info "  Best SDF level: $best_sdf"
    @info "  Best volume: $best_volume"  
    @info "  Best relative error: $(round(best_error * 100, digits=2))%"
    
    # Restore original positions and apply best solution
    for (idx, node_idx) in enumerate(surface_nodes)
        mesh.X[node_idx] = original_positions[idx]
        mesh.node_sdf[node_idx] = original_sdf_values[idx]
    end
    move_surface_nodes_to_sdf_level!(best_sdf)
    
    # Final volume assessment
    final_volume = calculate_current_mesh_volume()
    final_error = abs(final_volume - reference_volume) / reference_volume
    
    println("\nFinal Results:")
    println("  Reference volume: $reference_volume")
    println("  Final mesh volume: $final_volume") 
    println("  Final relative error: $(round(final_error * 100, digits=2))%")
    
    return final_error <= tolerance
end
