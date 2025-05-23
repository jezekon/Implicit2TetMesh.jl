"""
    slice_mesh_with_plane!(mesh::BlockMesh, plane::String, position::Float64, normal::Int; export_file::Union{String, Nothing}=nothing)

Slice a mesh with a specified plane, removing elements based on the normal direction.

Arguments:
- `mesh`: The BlockMesh object to be sliced
- `plane`: Cutting plane orientation, one of "x", "y", or "z"
- `position`: Relative position of the cut from 0.0 (beginning) to 1.0 (end)
- `normal`: Normal direction (1 = remove elements beyond plane, -1 = remove elements before plane)
- `export_file`: Optional filename to export the sliced mesh

Returns:
- The modified BlockMesh object with elements removed based on normal direction
"""
function slice_mesh_with_plane!(mesh::BlockMesh, plane::String, position::Float64, normal::Int; export_file::Union{String, Nothing}=nothing)
    # Validate input parameters
    if !(plane in ["z", "x", "y"])
        error("Invalid plane. Choose from: \"x\", \"y\", \"z\"")
    end
    
    if position < 0.0 || position > 1.0
        error("Position must be between 0 and 1")
    end
    
    if !(normal in [1, -1])
        error("Normal must be either 1 (remove beyond plane) or -1 (remove before plane)")
    end
    
    # Determine the axis perpendicular to the cutting plane
    cut_axis = if plane == "z" 
        3  # Z-axis is perpendicular to XY plane
    elseif plane == "x" 
        1  # X-axis is perpendicular to YZ plane
    else  # plane == "y"
        2  # Y-axis is perpendicular to XZ plane
    end
    
    # Calculate min and max coordinates along the cutting axis
    min_val = minimum(p[cut_axis] for p in mesh.X)
    max_val = maximum(p[cut_axis] for p in mesh.X)
    
    # Calculate the actual position of the cutting plane
    cut_pos = min_val + position * (max_val - min_val)
    
    # Determine removal direction based on normal
    direction_text = normal == 1 ? "beyond" : "before"
    @info "Cutting mesh with $(plane) plane at position $(position) (coordinate: $(cut_pos)), removing elements $(direction_text) the plane"
    
    # Count original elements and nodes for reporting
    orig_nodes = length(mesh.X)
    orig_elements = length(mesh.IEN)
    
    # Define a function to check if a point should be removed based on normal direction
    should_remove_point = if normal == 1
        p -> p[cut_axis] > cut_pos  # Remove points beyond the plane
    else  # normal == -1
        p -> p[cut_axis] < cut_pos  # Remove points before the plane
    end
    
    # Filter elements to keep only those with all nodes that should not be removed
    # This removes all tetrahedra that have at least one node in the removal zone
    mesh.IEN = [tet for tet in mesh.IEN if all(node_idx -> !should_remove_point(mesh.X[node_idx]), tet)]
    
    # Update mesh data structures to remove unused nodes and rebuild connectivity
    cleanup_unused_nodes!(mesh)  # Removes nodes that are no longer part of any element
    create_INE!(mesh)            # Rebuilds the inverse node-to-element connectivity
    
    # Count remaining elements and nodes
    remaining_nodes = length(mesh.X)
    remaining_elements = length(mesh.IEN)
    
    @info "Slice complete. Removed $(orig_nodes - remaining_nodes) nodes and $(orig_elements - remaining_elements) elements."
    
    # Export mesh if requested
    if export_file !== nothing
        export_mesh_vtu(mesh, export_file)
        @info "Exported sliced mesh to $(export_file)"
    end
    
    # return mesh
end

# Příklady použití:
# slice_mesh_with_plane!(mesh, "x", 0.6, 1, export_file="sliced_mesh_beyond.vtu")  # Smaže elementy za rovinou
# slice_mesh_with_plane!(mesh, "x", 0.6, -1, export_file="sliced_mesh_before.vtu") # Smaže elementy před rovinou
