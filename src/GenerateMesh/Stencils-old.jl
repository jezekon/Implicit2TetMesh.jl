# ----------------------------
# Helper function: Linear interpolation of the intersection with the zero level of SDF
# ----------------------------
function interpolate_zero(p1::SVector{3,Float64}, p2::SVector{3,Float64},
                          f1::Float64, f2::Float64, mesh::BlockMesh; tol=mesh.grid_tol, max_iter=20)::Int64
                          # coords of positive point, coords of negative point, sdf of positive point, sdf of negative point
                          
    pos1 = f1 >= -tol
    pos2 = f2 >= -tol
    # If both points have the same polarity according to tolerance, interpolation cannot be done correctly.
    if pos1 == pos2
        error("Both points have the same 'tolerance' polarity; one point must be close to zero (positive) and the other significantly negative.")
    end

    # Initialize interval: low and high - assuming p1 and p2 are ordered by f
    low, high = p1, p2
    f_low, f_high = f1, f2
    mid = low
    for iter in 1:max_iter
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
    sdf_of_iterp_point = eval_sdf(mesh, SVector{3, Float64}(p_key))
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

function apply_stencil(mesh::BlockMesh, tet::Vector{Int64})
  f = [ mesh.node_sdf[i] for i in tet ]
  np = count(x -> x > (0), f)
  if np > 0
    return [ tet ]
  else 
    return Vector{Vector{Int64}}()  # All negative - skip this tetrahedron
  end
end

# Function that traverses all tetrahedra in connectivity and applies complete stencils
function slice_ambiguous_tetrahedra!(mesh::BlockMesh)
  new_IEN = Vector{Vector{Int64}}()
  for tet in mesh.IEN
    # new_tets = apply_complete_stencil(mesh, tet)
    new_tets = apply_stencil(mesh, tet)
    for nt in new_tets
      t_sdf = [mesh.node_sdf[i] for i in nt]
      if any(x -> x >= 0, t_sdf)
        push!(new_IEN, nt)
      end
    end
  end
  mesh.IEN = new_IEN
  @info "After complete tetrahedral cutting: $(length(mesh.IEN)) tetrahedra"
end

"""
    remove_nodes_outside_isocontour!(mesh::BlockMesh, tol::Float64=mesh.grid_tol)

Removes all nodes with SDF values less than -tol and all tetrahedral elements 
that contain these nodes. This function helps clean up the mesh by removing
elements that lie outside the intended body defined by the zero isocontour.

# Arguments
- `mesh::BlockMesh`: The mesh to be cleaned
- `tol::Float64`: Tolerance value to determine which nodes are outside (default: mesh.grid_tol)
"""
function remove_nodes_outside_isocontour!(mesh::BlockMesh, tol::Float64=mesh.grid_tol*10000)
    @info "Removing nodes with SDF values less than -$tol and their connected elements..."
    
    # Identify nodes with SDF values less than -tol
    outside_nodes = Set{Int}()
    for (node_idx, sdf) in enumerate(mesh.node_sdf)
        if sdf < -tol
            push!(outside_nodes, node_idx)
        end
    end
    
    if isempty(outside_nodes)
        @info "No nodes found with SDF values less than -$tol."
        return mesh
    end
    
    @info "Found $(length(outside_nodes)) nodes with SDF values less than -$tol."
    
    # Identify elements containing these nodes using inverse connectivity
    # This is more efficient than checking each element individually
    outside_elements = Set{Int}()
    for node_idx in outside_nodes
        # Add all elements connected to this outside node
        union!(outside_elements, mesh.INE[node_idx])
    end
    
    # Count elements before removal
    original_element_count = length(mesh.IEN)
    
    # Keep only elements not marked for removal
    mesh.IEN = [element for (idx, element) in enumerate(mesh.IEN) if !(idx in outside_elements)]
    
    removed_count = original_element_count - length(mesh.IEN)
    @info "Removed $removed_count elements ($(round(removed_count/original_element_count*100, digits=2))%) containing nodes outside the isocontour."
    
    return mesh
end

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


# Function to warp nodes onto the isosurface based on element connectivity
function adjust_nodes_to_isosurface!(mesh::BlockMesh)
  @info "Starting node warping to isosurface..."
  
  # Track nodes that have been processed to avoid duplicating work
  processed_nodes = Set{Int}()
  nodes_moved = 0
  nodes_skipped = 0
  
  # Helper function to check if moving a node would cause element inversion
  function would_cause_inversion(node_idx, new_position)
      # Store original position
      original_position = mesh.X[node_idx]
      
      # Temporarily move the node
      mesh.X[node_idx] = new_position
      
      # Check all elements containing this node
      inversion_found = false
      for elem_idx in mesh.INE[node_idx]
          tet = mesh.IEN[elem_idx]
          vertices = [mesh.X[v] for v in tet]
          
          # Calculate Jacobian determinant (proportional to signed volume)
          jacobian_det = dot(
              cross(vertices[2] - vertices[1], vertices[3] - vertices[1]),
              vertices[4] - vertices[1]
          )
          
          if jacobian_det < 0
              inversion_found = true
              break
          end
      end
      
      # Restore original position
      mesh.X[node_idx] = original_position
      
      return inversion_found
  end
  
  # First pass: Find nodes with negative SDF and process them
  for node_idx in 1:length(mesh.X)
      # Skip already processed nodes
      node_idx in processed_nodes && continue
      
      # Only process nodes with negative SDF values
      if mesh.node_sdf[node_idx] < 0
          # Get all elements containing this node
          elements = mesh.INE[node_idx]
          
          # Skip nodes that aren't part of any elements
          isempty(elements) && continue
          
          # Try to find an element with exactly 3 negative nodes
          found_3neg_element = false
          edge_to_warp = nothing
          
          for elem_idx in elements
              tet = mesh.IEN[elem_idx]
              # Count negative nodes in this tetrahedron
              neg_nodes = filter(i -> mesh.node_sdf[i] < 0, tet)
              
              if length(neg_nodes) == 3 && node_idx in neg_nodes
                  # This element has exactly 3 negative nodes including our target node
                  found_3neg_element = true
                  
                  # Find the positive node
                  pos_node = first(filter(i -> mesh.node_sdf[i] >= 0, tet))
                  
                  # Calculate the SDF values
                  pos_sdf = mesh.node_sdf[pos_node]
                  neg_sdf = mesh.node_sdf[node_idx]
                  
                  # Check if the edge between positive and negative node crosses zero
                  if sign(pos_sdf) != sign(neg_sdf)
                      # Calculate approximate zero crossing using linear interpolation
                      t = pos_sdf / (pos_sdf - neg_sdf)
                      zero_point = mesh.X[pos_node] + t * (mesh.X[node_idx] - mesh.X[pos_node])
                      
                      # Check if moving to this point would cause inversion
                      if !would_cause_inversion(node_idx, zero_point)
                          # Found a suitable edge to warp along
                          edge_to_warp = (pos_node, node_idx)
                          break
                      end
                  end
              end
          end
          
          # If no suitable edge found yet, try all edges connecting to positive nodes
          if edge_to_warp === nothing
              for elem_idx in elements
                  tet = mesh.IEN[elem_idx]
      
                  # Look for edges connecting our negative node to positive nodes
                  for other_node in tet
                      # Skip if it's the same node
                      other_node == node_idx && continue
          
                      # Check if the other node has positive SDF (crosses zero level)
                      if mesh.node_sdf[other_node] > 0
                          # Calculate approximate zero crossing
                          pos_sdf = mesh.node_sdf[other_node]
                          neg_sdf = mesh.node_sdf[node_idx]
                          
                          # Make sure signs are opposite (to cross the zero level)
                          if sign(pos_sdf) == sign(neg_sdf)
                              continue
                          end
                          
                          t = pos_sdf / (pos_sdf - neg_sdf)
                          zero_point = mesh.X[other_node] + t * (mesh.X[node_idx] - mesh.X[other_node])
                          
                          # Check if moving to this point would cause inversion
                          if !would_cause_inversion(node_idx, zero_point)
                              # Found a suitable edge to warp along
                              edge_to_warp = (other_node, node_idx)
                              break
                          end
                      end
                  end
      
                  # If we found a suitable edge, no need to check more elements
                  edge_to_warp !== nothing && break
              end
          end
          
          # If we found an edge to warp along, perform the warping
          if edge_to_warp !== nothing
              pos_node, neg_node = edge_to_warp
              
              # Interpolate to find the zero crossing
              try
                  new_node_idx = interpolate_zero(
                      mesh.X[pos_node],  # coordinates of positive vertex
                      mesh.X[neg_node],  # coordinates of negative vertex
                      mesh.node_sdf[pos_node],  # SDF value of positive vertex
                      mesh.node_sdf[neg_node],  # SDF value of negative vertex
                      mesh
                  )
                  
                  # Update the node's position and SDF value if a new node wasn't created
                  if new_node_idx == neg_node
                      mesh.node_sdf[neg_node] = 0.0
                      nodes_moved += 1
                  else
                      # If a new node was created, update all elements that contain the old node
                      for elem_idx in elements
                          tet = mesh.IEN[elem_idx]
                          for i in 1:length(tet)
                              if tet[i] == neg_node
                                  mesh.IEN[elem_idx][i] = new_node_idx
                              end
                          end
                      end
                      nodes_moved += 1
                  end
              catch e
                  @warn "Error during interpolation for node $node_idx: $e"
                  nodes_skipped += 1
              end
          else
              # No suitable edge was found - all would cause inversion
              @warn "Node $node_idx couldn't be moved to isosurface - all potential movements would cause element inversion"
              nodes_skipped += 1
          end
          
          # Mark this node as processed
          push!(processed_nodes, node_idx)
      end
  end
  
  @info "Node warping complete. $nodes_moved nodes were moved to the isosurface. $nodes_skipped nodes were skipped to prevent element inversion."
 
  # Remove nodes outside the isocontour and elements that contain them
  # remove_nodes_outside_isocontour!(mesh)

  # Remove any inverted elements (with negative Jacobian determinant)
  remove_inverted_elements!(mesh)

  # After warping, we need to update the mesh connectivity
  update_connectivity!(mesh)
end


