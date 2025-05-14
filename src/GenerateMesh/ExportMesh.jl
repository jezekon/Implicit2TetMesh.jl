# ----------------------------
# Export mesh to VTU (unchanged logic, minor type adjustments)
# ----------------------------
function export_mesh_vtu(mesh::BlockMesh, filename::String)
  npoints = length(mesh.X)
  points = zeros(Float64, 3, npoints)
  @inbounds for i in 1:npoints
    points[:, i] = mesh.X[i]
  end
  cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
  vtkfile = vtk_grid(filename, points, cells)
  vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
  vtk_save(vtkfile)
end

# ----------------------------
# Enhanced VTU export with element volume and quality metrics
# ----------------------------
function export_mesh_vtu_quality(mesh::BlockMesh, filename::String)
    # Prepare points data
    npoints = length(mesh.X)
    points = zeros(Float64, 3, npoints)
    @inbounds for i in 1:npoints
        points[:, i] = mesh.X[i]
    end
    
    # Create cells
    cells = [MeshCell(VTKCellTypes.VTK_TETRA, tet) for tet in mesh.IEN]
    
    # Create VTK grid
    vtkfile = vtk_grid(filename, points, cells)
    
    # Add point data (SDF values)
    vtk_point_data(vtkfile, mesh.node_sdf, "sdf")
    
    # Calculate element volumes
    nelements = length(mesh.IEN)
    volumes = zeros(Float64, nelements)
    quality = zeros(Float64, nelements)
    
    # Calculate volume and quality for each tetrahedron
    @inbounds for i in 1:nelements
        vertices = [mesh.X[node_idx] for node_idx in mesh.IEN[i]]
        
        # Calculate volume using determinant method
        a = vertices[2] - vertices[1]
        b = vertices[3] - vertices[1]
        c = vertices[4] - vertices[1]
        
        # Volume = (1/6) * |a · (b × c)|
        volumes[i] = abs(dot(a, cross(b, c))) / 6.0
        
        # Calculate quality using mean ratio metric
        # This is a scale-invariant measure from 0 (degenerate) to 1 (perfect)
        
        # Calculate all edge vectors
        edges = [
            vertices[2] - vertices[1], vertices[3] - vertices[1], vertices[4] - vertices[1],
            vertices[3] - vertices[2], vertices[4] - vertices[2],
            vertices[4] - vertices[3]
        ]
        
        # Calculate sum of squared edge lengths
        sum_squared_lengths = sum(dot(e, e) for e in edges)
        
        # Mean ratio quality metric: 12.0 * (3.0 * volume)^(2.0/3.0) / sum(edge_lengths^2)
        # This evaluates to 1.0 for a perfect equilateral tetrahedron
        if volumes[i] > 1e-10  # Avoid division by zero and numerical issues
            quality[i] = 12.0 * (3.0 * volumes[i])^(2.0/3.0) / sum_squared_lengths
        else
            quality[i] = 0.0  # Degenerate element
        end
    end
    
    # Add cell data (volume and quality metrics)
    vtk_cell_data(vtkfile, volumes, "volume")
    vtk_cell_data(vtkfile, quality, "quality")
    
    # Save the VTK file
    vtk_save(vtkfile)
    
    # Return summary statistics
    return Dict(
        "min_volume" => minimum(volumes),
        "max_volume" => maximum(volumes),
        "avg_volume" => sum(volumes) / nelements,
        "min_quality" => minimum(quality),
        "max_quality" => maximum(quality),
        "avg_quality" => sum(quality) / nelements
    )
end
