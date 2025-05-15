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
# export_mesh_vtu_quality


# ----------------------------
# Enhanced VTU export with element volume, quality and Jacobian determinant
# Optimized for easy visualization of problematic elements in Paraview
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
    
    # Calculate element volumes, qualities, and Jacobian determinants
    nelements = length(mesh.IEN)
    volumes = zeros(Float64, nelements)
    qualities = zeros(Float64, nelements)
    jacobian_dets = zeros(Float64, nelements)  # Raw Jacobian determinants
    jacobian_quality = zeros(Float64, nelements)  # Normalized quality based on Jacobian
    volume_quality = zeros(Float64, nelements)  # Normalized quality based on volume
    problematic_elements = zeros(Int, nelements)  # Indicator for problematic elements
    
    # Ideal dihedral angle for a perfect tetrahedron (approximately arccos(-1/3))
    ideal_angle = 70.53
    
    # Calculate volume, quality, and Jacobian determinant for each tetrahedron
    @inbounds for i in 1:nelements
        # Get vertices of the tetrahedron
        vertices = [mesh.X[node_idx] for node_idx in mesh.IEN[i]]
        
        # Calculate edge vectors
        a = vertices[2] - vertices[1]
        b = vertices[3] - vertices[1]
        c = vertices[4] - vertices[1]
        
        # Calculate Jacobian determinant (proportional to signed volume)
        det_value = dot(a, cross(b, c))
        jacobian_dets[i] = det_value  # Store the raw determinant
        
        # Calculate volume (absolute value of 1/6 of the determinant)
        volumes[i] = abs(det_value) / 6.0
        
        # Calculate face normals (outward facing)
        normals = [
            normalize(cross(vertices[2] - vertices[1], vertices[3] - vertices[1])),  # Face 123
            normalize(cross(vertices[2] - vertices[1], vertices[4] - vertices[1])),  # Face 124
            normalize(cross(vertices[3] - vertices[1], vertices[4] - vertices[1])),  # Face 134
            normalize(cross(vertices[3] - vertices[2], vertices[4] - vertices[2]))   # Face 234
        ]
        
        # Ensure normals are consistently outward-facing
        face_centers = [
            (vertices[1] + vertices[2] + vertices[3]) / 3.0,  # Center of face 123
            (vertices[1] + vertices[2] + vertices[4]) / 3.0,  # Center of face 124
            (vertices[1] + vertices[3] + vertices[4]) / 3.0,  # Center of face 134
            (vertices[2] + vertices[3] + vertices[4]) / 3.0   # Center of face 234
        ]
        
        opposite_vertices = [vertices[4], vertices[3], vertices[2], vertices[1]]
        
        # Flip normals if needed to ensure they point outward
        for j in 1:4
            if dot(normals[j], opposite_vertices[j] - face_centers[j]) > 0
                normals[j] = -normals[j]
            end
        end
        
        # Calculate all 6 dihedral angles (in degrees)
        dihedral_angles = [
            acos(-clamp(dot(normals[i], normals[j]), -1.0, 1.0)) * 180 / π
            for i in 1:3 for j in i+1:4
        ]
        
        # Find the angle that deviates most from the ideal angle
        deviations = abs.(dihedral_angles .- ideal_angle)
        max_deviation = maximum(deviations)
        
        # Calculate quality based on maximum deviation from ideal angle
        # Quality ranges from 0 (degenerate) to 1 (perfect)
        if max_deviation >= ideal_angle  # Angle close to 0° or 180°
            qualities[i] = 0.0
        else
            # Linear mapping from maximum deviation to quality
            # When deviation = 0, quality = 1
            # When deviation = ideal_angle, quality = 0
            qualities[i] = 1.0 - (max_deviation / ideal_angle)
        end
    end
    
    # Calculate threshold values for marking problematic elements
    # For Jacobian determinant: Negative values are critically bad, small positive values are concerning
    # Find the median positive Jacobian determinant as a reference
    positive_dets = filter(d -> d > 0, jacobian_dets)
    if !isempty(positive_dets)
        median_positive_det = median(positive_dets)
        threshold_det = median_positive_det * 0.01  # 1% of median is our threshold for concerning
    else
        threshold_det = 1.0  # Fallback if no positive determinants (unlikely)
    end
    
    # For volume: Very small volumes relative to the median are concerning
    median_volume = median(volumes)
    threshold_volume = median_volume * 0.01  # 1% of median is our threshold
    
    # Process each element to compute quality indicators
    @inbounds for i in 1:nelements
        # Jacobian quality: 0 for negative, between 0-1 for small positive, 1 for good
        if jacobian_dets[i] < 0
            jacobian_quality[i] = 0.0  # Critical issue - negative determinant
            problematic_elements[i] = 1  # Mark as problematic
        elseif jacobian_dets[i] < threshold_det
            # Scale quality between 0-0.5 for small positive determinants
            # (0 = 0, threshold = 0.5)
            jacobian_quality[i] = 0.5 * (jacobian_dets[i] / threshold_det)
            problematic_elements[i] = 1  # Mark as problematic
        else
            # Scale quality between 0.5-1 for good determinants
            # (threshold = 0.5, max = 1)
            jacobian_quality[i] = 0.5 + 0.5 * min(jacobian_dets[i] / (10 * threshold_det), 1.0)
        end
        
        # Volume quality: 0 for zero volume, scales up to 1 for good volumes
        if volumes[i] < threshold_volume
            # Scale quality between 0-0.5 for small volumes
            volume_quality[i] = 0.5 * (volumes[i] / threshold_volume)
            problematic_elements[i] = 1  # Mark as problematic
        else
            # Scale quality between 0.5-1 for good volumes
            volume_quality[i] = 0.5 + 0.5 * min(volumes[i] / (10 * threshold_volume), 1.0)
        end
        
        # Also check if quality based on dihedral angles is very poor
        if qualities[i] < 0.2  # Quality less than 0.2 is concerning
            problematic_elements[i] = 1  # Mark as problematic
        end
    end
    
    # Find element with worst quality
    worst_quality_idx = argmin(qualities)
    worst_quality = qualities[worst_quality_idx]
    
    # Check for negative Jacobian determinants (inverted elements)
    negative_jacobian_count = count(j -> j < 0, jacobian_dets)
    if negative_jacobian_count > 0
        println("WARNING: Found $negative_jacobian_count elements with negative Jacobian determinant!")
        # Find element with most negative Jacobian
        worst_jacobian_idx = argmin(jacobian_dets)
        println("  Most inverted element index: ", worst_jacobian_idx)
        println("  Jacobian determinant value: ", jacobian_dets[worst_jacobian_idx])
        println("  Element nodes: ", mesh.IEN[worst_jacobian_idx])
    else
        println("All elements have positive Jacobian determinants (valid for FEM).")
    end
    
    # Count problematic elements
    problem_count = sum(problematic_elements)
    println("Identified $problem_count problematic elements out of $nelements total ($(round(100*problem_count/nelements, digits=2))%).")
    
    # Print information about element with worst quality
    println("Element with worst quality:")
    println("  Element index: ", worst_quality_idx)
    println("  Quality value: ", worst_quality)
    println("  Element volume: ", volumes[worst_quality_idx])
    println("  Jacobian determinant: ", jacobian_dets[worst_quality_idx])
    println("  Element nodes: ", mesh.IEN[worst_quality_idx])
    
    # Add cell data - first the raw data
    vtk_cell_data(vtkfile, volumes, "volume")
    vtk_cell_data(vtkfile, jacobian_dets, "jacobian_det")
    vtk_cell_data(vtkfile, qualities, "angle_quality")
    
    # Add the normalized quality indicators designed for visualization
    vtk_cell_data(vtkfile, jacobian_quality, "jacobian_quality")
    vtk_cell_data(vtkfile, volume_quality, "volume_quality")
    vtk_cell_data(vtkfile, problematic_elements, "problematic")
    
    # Save the VTK file
    vtk_save(vtkfile)
    
    # Print usage tips for Paraview
    println("\nParaview visualization tips:")
    println("1. Load the file and select 'Surface With Edges' for better element visibility")
    println("2. Use 'problematic' field with threshold filter set to 1 to isolate problematic elements")
    println("3. For Jacobian quality, use a diverging color map with 0.5 at center:")
    println("   - Red (0.0): Critical - negative Jacobian determinant")
    println("   - Yellow (0.5): Borderline acceptable")
    println("   - Green (1.0): Good quality")
    println("4. Same approach for volume_quality visualization")
    
    # Return summary statistics
    return Dict(
        "min_volume" => minimum(volumes),
        "max_volume" => maximum(volumes),
        "avg_volume" => sum(volumes) / nelements,
        "min_quality" => minimum(qualities),
        "max_quality" => maximum(qualities),
        "avg_quality" => sum(qualities) / nelements,
        "min_jacobian_det" => minimum(jacobian_dets),
        "max_jacobian_det" => maximum(jacobian_dets),
        "avg_jacobian_det" => sum(jacobian_dets) / nelements,
        "negative_jacobian_count" => negative_jacobian_count,
        "problematic_element_count" => problem_count,
        "worst_quality_element" => worst_quality_idx
    )
end
