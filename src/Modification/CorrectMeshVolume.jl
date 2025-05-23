"""
    detect_boundary_nodes(mesh::BlockMesh) -> Set{Int}

Detekuje uzly na hranici sítě (uzly patřící k vnějším stěnám).
"""
function detect_boundary_nodes(mesh::BlockMesh)
    # Slovník pro počítání, kolikrát se jednotlivé trojúhelníkové stěny vyskytují
    face_count = Dict{Tuple{Int,Int,Int}, Int}()
    
    # Procházej všechny tetrahedrální elementy a jejich stěny
    for tet in mesh.IEN
        # Každý tetrahedr má 4 stěny
        faces = [
            tuple(sort([tet[1], tet[2], tet[3]])...),  # Stěna 1-2-3
            tuple(sort([tet[1], tet[2], tet[4]])...),  # Stěna 1-2-4
            tuple(sort([tet[1], tet[3], tet[4]])...),  # Stěna 1-3-4
            tuple(sort([tet[2], tet[3], tet[4]])...)   # Stěna 2-3-4
        ]
        
        for face in faces
            face_count[face] = get(face_count, face, 0) + 1
        end
    end
    
    # Hranični stěny se vyskytují pouze jednou
    boundary_faces = [face for (face, count) in face_count if count == 1]
    
    # Sesbírej všechny uzly z hraničních stěn
    boundary_nodes = Set{Int}()
    for face in boundary_faces
        union!(boundary_nodes, face)
    end
    
    return boundary_nodes
end

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
    correct_mesh_volume!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3}; 
                        tolerance::Float64=0.0005, max_iterations::Int=50)

Koriguje objem tetrahedrální sítě tak, aby přesně odpovídal objemu implicitní geometrie.
Používá metodu půlení intervalu pro iterativní přiblížení k cílovému objemu.
"""
function correct_mesh_volume!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3}; 
                             tolerance::Float64=0.0005, max_iterations::Int=50)
    @info "Correcting mesh volume to match SDF geometry..."
    
    # 1. Vypočítej referenční objem z SDF dat
    reference_volume = Float64(calculate_volume_from_sdf(fine_sdf, fine_grid))
    @info "Reference volume from SDF: $reference_volume"
    
    # 2. Detekuj vnější uzly sítě
    boundary_nodes = detect_boundary_nodes(mesh)
    @info "Detected $(length(boundary_nodes)) boundary nodes"
    
    # 3. Uložení původních pozic hraničních uzlů
    original_positions = Dict{Int, SVector{3,Float64}}()
    for node_idx in boundary_nodes
        original_positions[node_idx] = mesh.X[node_idx]
    end
    
    # 4. Inicializuj parametry pro půlení intervalu
    initial_step = 0.1  # Zvětšený počáteční krok
    lower_bound = -initial_step
    upper_bound = initial_step
    
    # Výpočet počátečního objemu sítě
    mesh_volume = calculate_mesh_volume(mesh)
    @info "Initial mesh volume: $mesh_volume"
    
    relative_error = abs(mesh_volume - reference_volume) / reference_volume
    @info "Initial relative error: $(relative_error * 100)%"
    
    if relative_error <= tolerance
        @info "Mesh volume already within tolerance"
        return true
    end
    
    # Inicializace pro půlení intervalu
    best_level = 0.0
    best_error = relative_error
    
    # 5. Hlavní iterační smyčka s půlením intervalu
    for iteration in 1:max_iterations
        # Určení nové úrovne SDF pro test
        current_level = (lower_bound + upper_bound) / 2.0
        
        @info "Iteration $iteration: Testing SDF level $current_level"
        
        # Resetuj pozice uzlů na původní
        restore_boundary_nodes!(mesh, boundary_nodes, original_positions)
        
        # Posuň hranční uzly na novou úroveň SDF
        successful_moves = move_boundary_nodes_to_sdf_level!(mesh, boundary_nodes, current_level)
        
        if successful_moves == 0
            @warn "No nodes could be moved to SDF level $current_level"
            # Zmenši interval
            if mesh_volume < reference_volume
                upper_bound = current_level
            else
                lower_bound = current_level
            end
            continue
        end
        
        # Vypočítaj nový objem sítě
        new_mesh_volume = calculate_mesh_volume(mesh)
        new_relative_error = abs(new_mesh_volume - reference_volume) / reference_volume
        
        @info "  New mesh volume: $new_mesh_volume"
        @info "  Relative error: $(new_relative_error * 100)%"
        @info "  Successfully moved: $successful_moves nodes"
        
        # Aktualizuj nejlepší řešení
        if new_relative_error < best_error
            best_error = new_relative_error
            best_level = current_level
        end
        
        # Kontrola konvergence
        if new_relative_error <= tolerance
            @info "Volume correction converged after $iteration iterations"
            @info "Final volume: $new_mesh_volume (target: $reference_volume)"
            return true
        end
        
        # Logika půlení intervalu - OPRAVENÁ
        if new_mesh_volume < reference_volume
            # Objem sítě je menší → potřebujeme větší objem → posunout uzly směrem ven (pozitivní SDF)
            lower_bound = current_level
        else
            # Objem sítě je větší → potřebujeme menší objem → posunout uzly směrem dovnitř (negativní SDF)
            upper_bound = current_level
        end
        
        # Kontrola, zda je interval příliš malý pro další iterace
        if abs(upper_bound - lower_bound) < 1e-8
            @warn "Interval became too small, stopping iterations"
            break
        end
        
        mesh_volume = new_mesh_volume
    end
    
    # Použij nejlepší nalezené řešení
    @info "Applying best solution with SDF level $best_level (error: $(best_error * 100)%)"
    restore_boundary_nodes!(mesh, boundary_nodes, original_positions)
    move_boundary_nodes_to_sdf_level!(mesh, boundary_nodes, best_level)
    
    @warn "Volume correction did not converge within $max_iterations iterations"
    @info "Final relative error: $(best_error * 100)%"
    return false
end

"""
    restore_boundary_nodes!(mesh::BlockMesh, boundary_nodes::Set{Int}, original_positions::Dict{Int, SVector{3,Float64}})

Obnoví původní pozice hraničních uzlů.
"""
function restore_boundary_nodes!(mesh::BlockMesh, boundary_nodes::Set{Int}, original_positions::Dict{Int, SVector{3,Float64}})
    for node_idx in boundary_nodes
        mesh.X[node_idx] = original_positions[node_idx]
        mesh.node_sdf[node_idx] = eval_sdf(mesh, mesh.X[node_idx])
    end
end

"""
    move_boundary_nodes_to_sdf_level!(mesh::BlockMesh, boundary_nodes::Set{Int}, target_sdf::Float64; 
                                     max_iterations::Int=10, tolerance::Float64=1e-4) -> Int

Posune hranční uzly na specifikovanou úroveň SDF. Vrací počet úspěšně posunutých uzlů.
"""
function move_boundary_nodes_to_sdf_level!(mesh::BlockMesh, boundary_nodes::Set{Int}, target_sdf::Float64; 
                                          max_iterations::Int=10, tolerance::Float64=1e-4)
    nodes_moved = 0
    
    for node_idx in boundary_nodes
        success = move_node_to_sdf_level!(mesh, node_idx, target_sdf, max_iterations, tolerance)
        if success
            nodes_moved += 1
        end
    end
    
    return nodes_moved
end

"""
    move_node_to_sdf_level!(mesh::BlockMesh, node_idx::Int, target_sdf::Float64, 
                           max_iterations::Int, tolerance::Float64) -> Bool

Posune jednotlivý uzel na specifikovanou úroveň SDF s omezením vzdálenosti.
"""
function move_node_to_sdf_level!(mesh::BlockMesh, node_idx::Int, target_sdf::Float64, 
                                max_iterations::Int, tolerance::Float64)
    original_position = mesh.X[node_idx]
    current_position = original_position
    max_displacement = mesh.grid_step * 0.5  # Omezení maximálního posuvu
    
    for iteration in 1:max_iterations
        # Vypočítaj aktuální SDF hodnotu
        current_sdf = eval_sdf(mesh, current_position)
        
        # Kontrola konvergence
        if abs(current_sdf - target_sdf) < tolerance
            # Kontrola, zda posuv není příliš velký
            displacement = norm(current_position - original_position)
            if displacement <= max_displacement
                mesh.X[node_idx] = current_position
                mesh.node_sdf[node_idx] = current_sdf
                return true
            else
                # Omezení posuvu
                direction = normalize(current_position - original_position)
                mesh.X[node_idx] = original_position + max_displacement * direction
                mesh.node_sdf[node_idx] = eval_sdf(mesh, mesh.X[node_idx])
                return true
            end
        end
        
        # Vypočítaj gradient SDF
        gradient = compute_gradient(mesh, current_position)
        gradient_norm_sq = sum(abs2, gradient)
        
        # Kontrola, zda gradient není příliš malý
        if gradient_norm_sq < 1e-12
            return false
        end
        
        # Newton-Raphsonův krok s adaptivním krokem
        step_size = (current_sdf - target_sdf) / gradient_norm_sq
        step_size = clamp(step_size, -mesh.grid_step * 0.1, mesh.grid_step * 0.1)  # Omezení kroku
        
        new_position = current_position - step_size * gradient
        
        # Kontrola, zda se nepřesáhnula maximální vzdálenost
        total_displacement = norm(new_position - original_position)
        if total_displacement > max_displacement
            direction = normalize(new_position - original_position)
            new_position = original_position + max_displacement * direction
        end
        
        current_position = new_position
    end
    
    # Pokud se nepodařilo konvergovat, použij alespoň částečný posuv
    displacement = norm(current_position - original_position)
    if displacement <= max_displacement
        mesh.X[node_idx] = current_position
        mesh.node_sdf[node_idx] = eval_sdf(mesh, mesh.X[node_idx])
        return true
    end
    
    return false
end

"""
    adaptive_volume_correction!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3}; 
                               tolerance::Float64=0.0005, max_iterations::Int=30)

Alternativní přístup používající adaptivní metodu pro korekci objemu.
"""
function adaptive_volume_correction!(mesh::BlockMesh, fine_sdf::Array{Float32,3}, fine_grid::Array{Vector{Float32},3}; 
                                   tolerance::Float64=0.000005, max_iterations::Int=30)
    @info "Starting adaptive volume correction..."
    
    reference_volume = Float64(calculate_volume_from_sdf(fine_sdf, fine_grid))
    boundary_nodes = detect_boundary_nodes(mesh)
    
    # Uložení původních pozic
    original_positions = Dict{Int, SVector{3,Float64}}()
    for node_idx in boundary_nodes
        original_positions[node_idx] = mesh.X[node_idx]
    end
    
    current_volume = calculate_mesh_volume(mesh)
    relative_error = abs(current_volume - reference_volume) / reference_volume
    
    @info "Initial state: mesh_volume=$current_volume, reference_volume=$reference_volume, error=$(relative_error*100)%"
    
    if relative_error <= tolerance
        return true
    end
    
    # Adaptivní přístup
    step_size = 0.02
    direction = current_volume < reference_volume ? 1.0 : -1.0  # 1 = ven, -1 = dovnitř
    
    best_volume = current_volume
    best_error = relative_error
    best_step = 0.0
    
    for iteration in 1:max_iterations
        test_sdf_level = direction * step_size * iteration
        
        # Resetuj pozice
        restore_boundary_nodes!(mesh, boundary_nodes, original_positions)
        
        # Aplikuj korekci
        successful_moves = move_boundary_nodes_to_sdf_level!(mesh, boundary_nodes, test_sdf_level)
        
        if successful_moves > 0
            new_volume = calculate_mesh_volume(mesh)
            new_error = abs(new_volume - reference_volume) / reference_volume
            
            @info "Iteration $iteration: SDF_level=$test_sdf_level, volume=$new_volume, error=$(new_error*100)%, moved=$successful_moves"
            
            if new_error < best_error
                best_error = new_error
                best_volume = new_volume
                best_step = test_sdf_level
            end
            
            if new_error <= tolerance
                @info "Adaptive correction converged!"
                return true
            end
        end
    end
    
    # Aplikuj nejlepší řešení
    restore_boundary_nodes!(mesh, boundary_nodes, original_positions)
    move_boundary_nodes_to_sdf_level!(mesh, boundary_nodes, best_step)
    
    @info "Best solution: SDF_level=$best_step, error=$(best_error*100)%"
    return best_error <= tolerance
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
