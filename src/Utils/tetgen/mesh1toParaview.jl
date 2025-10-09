"""
Convert surface (stl) mesh to volume (tetrahedra) mesh using Tetgen:
tetgen -p model.stl

This scipt converts mesh1 data to vtu.

"""

using WriteVTK

function tetgen_to_vtk(base_filename, output_filename)
    # Načtení uzlů (.node soubor)
    node_file = open("$base_filename.node", "r")
    lines = readlines(node_file)
    close(node_file)

    # Zpracování hlavičky
    header = split(lines[1])
    num_points = parse(Int, header[1])

    # Alokace pole pro body - potřebujeme znát přesně počet bodů
    points = zeros(Float64, 3, num_points)

    # Zjištění, zda indexování začíná od 0 nebo od 1 podle prvního bodu
    first_point_line = split(lines[2])
    first_index = parse(Int, first_point_line[1])
    index_offset = (first_index == 0) ? 1 : 0  # Pokud indexy začínají od 0, přidáme offset 1

    # Načtení bodů
    for i = 1:num_points
        values = split(lines[i+1])
        # Převedení indexu z TetGen na index do pole (0-based -> 1-based pokud je potřeba)
        point_idx = parse(Int, values[1]) + index_offset
        if point_idx <= num_points
            points[1, point_idx] = parse(Float64, values[2])
            points[2, point_idx] = parse(Float64, values[3])
            points[3, point_idx] = parse(Float64, values[4])
        else
            println("Varování: Index bodu $point_idx je mimo rozsah (max $num_points)")
        end
    end

    # Načtení prvků (.ele soubor)
    ele_file = open("$base_filename.ele", "r")
    lines = readlines(ele_file)
    close(ele_file)

    # Zpracování hlavičky
    header = split(lines[1])
    num_tets = parse(Int, header[1])
    nodes_per_tet = length(split(lines[2])) - 1  # Počet uzlů na jeden tetrahedron

    # Načtení tetraedrů
    cells = Array{MeshCell}(undef, num_tets)
    for i = 1:num_tets
        values = split(lines[i+1])
        # Přeskočení indexu prvku, načtení indexů vrcholů tetraedru
        tet_indices = [parse(Int, values[j]) + index_offset for j = 2:(nodes_per_tet+1)]
        cells[i] = MeshCell(VTKCellTypes.VTK_TETRA, tet_indices)
    end

    # Vytvoření VTK souboru
    vtk_grid(output_filename, points, cells) do vtk
        # Zde můžete přidat další atributy, pokud jsou dostupné
    end

    println("Převod dokončen. Soubor uložen jako $(output_filename).vtu")
end

# Použití pro vaše konkrétní soubory
base_filename = "cantilever_beam_interp_smooth1_iso-tri_mesh.1"
output_filename = "cantilever_beam_tet_interp"
tetgen_to_vtk(base_filename, output_filename)
