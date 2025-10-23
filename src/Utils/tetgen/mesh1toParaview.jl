"""
Convert surface (STL) mesh to volume (tetrahedra) mesh using Tetgen:
tetgen -p model.stl

This script converts mesh1 data to VTU format.
"""

using WriteVTK

function tetgen_to_vtk(base_filename, output_filename)
    # Load nodes from .node file
    node_file = open("$base_filename.node", "r")
    lines = readlines(node_file)
    close(node_file)

    # Parse header
    header = split(lines[1])
    num_points = parse(Int, header[1])

    # Allocate array for points - we need to know exact number of points
    points = zeros(Float64, 3, num_points)

    # Determine if indexing starts from 0 or 1 based on first point
    first_point_line = split(lines[2])
    first_index = parse(Int, first_point_line[1])
    index_offset = (first_index == 0) ? 1 : 0  # If indices start from 0, add offset of 1

    # Load points
    for i = 1:num_points
        values = split(lines[i+1])
        # Convert index from TetGen to array index (0-based -> 1-based if needed)
        point_idx = parse(Int, values[1]) + index_offset
        if point_idx <= num_points
            points[1, point_idx] = parse(Float64, values[2])
            points[2, point_idx] = parse(Float64, values[3])
            points[3, point_idx] = parse(Float64, values[4])
        else
            println("Warning: Point index $point_idx is out of range (max $num_points)")
        end
    end

    # Load elements from .ele file
    ele_file = open("$base_filename.ele", "r")
    lines = readlines(ele_file)
    close(ele_file)

    # Parse header
    header = split(lines[1])
    num_tets = parse(Int, header[1])
    nodes_per_tet = length(split(lines[2])) - 1  # Number of nodes per tetrahedron

    # Load tetrahedra
    cells = Array{MeshCell}(undef, num_tets)
    for i = 1:num_tets
        values = split(lines[i+1])
        # Skip element index, read vertex indices of tetrahedron
        tet_indices = [parse(Int, values[j]) + index_offset for j = 2:(nodes_per_tet+1)]
        cells[i] = MeshCell(VTKCellTypes.VTK_TETRA, tet_indices)
    end

    # Create VTK file
    vtk_grid(output_filename, points, cells) do vtk
        # Additional attributes can be added here if available
    end

    println("Conversion complete. File saved as $(output_filename).vtu")
end

# Usage for your specific files
base_filename = "cantilever_beam_interp_smooth1_iso-tri_mesh.1"
output_filename = "cantilever_beam_tet_interp"
tetgen_to_vtk(base_filename, output_filename)
