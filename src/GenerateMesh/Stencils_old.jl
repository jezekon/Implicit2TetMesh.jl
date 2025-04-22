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
# function slice_ambiguous_tetrahedra!(mesh::BlockMesh)
#   new_IEN = Vector{Vector{Int64}}()
#   for tet in mesh.IEN
#     # new_tets = apply_complete_stencil(mesh, tet)
#     new_tets = apply_stencil(mesh, tet)
#     for nt in new_tets
#       t_sdf = [mesh.node_sdf[i] for i in nt]
#       if any(x -> x >= 0, t_sdf)
#         push!(new_IEN, nt)
#       end
#     end
#   end
#   mesh.IEN = new_IEN
#   @info "After complete tetrahedral cutting: $(length(mesh.IEN)) tetrahedra"
# end

function slice_ambiguous_tetrahedra!(mesh::BlockMesh)
    @info "Slicing tetrahedra using trim_spikes logic..."
    new_IEN = Vector{Vector{Int64}}()
    sizehint!(new_IEN, length(mesh.IEN)) # Pre-allocate estimate

    # Keep track of edges cut during this process
    cut_map = Dict{Tuple{Int, Int}, Int}()

    original_node_count = length(mesh.X)
    original_tet_count = length(mesh.IEN)

    # Iterate through a copy of IEN if modifying mesh.X/node_sdf inside loop
    # Or ensure cut_edge! handles node addition carefully.
    # Let's iterate through original IEN and build new_IEN.
    current_IEN = mesh.IEN
    mesh.IEN = Vector{Vector{Int64}}() # Clear current IEN temporarily

    for tet in current_IEN
        # Apply the trim_spikes stencil
        # Pass mesh.node_sdf explicitly to cut_edge! via apply_stencil_trim_spikes!
        resulting_tets = apply_stencil_trim_spikes!(mesh, tet, cut_map)

        # Add the resulting valid tetrahedra to the new list
        for nt in resulting_tets
            # Optional: Add a check here to ensure the tet is valid (e.g., positive volume)
            # before adding it to new_IEN. The C++ code does this later.
            push!(new_IEN, nt)
        end
    end

    mesh.IEN = new_IEN # Assign the newly generated tetrahedra list
    new_node_count = length(mesh.X)
    new_tet_count = length(mesh.IEN)

    @info "After trim_spikes slicing: $(new_tet_count) tetrahedra (added $(new_tet_count - original_tet_count)), $(new_node_count) nodes (added $(new_node_count - original_node_count))"

    # Important: After adding nodes and changing connectivity,
    # the mesh needs cleanup and INE update.
    # update_connectivity!(mesh) # Call this *after* slicing is complete
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

# Stencils.jl
# (Zachovejte existující includes a definice struktur)

# using StaticArrays
# using LinearAlgebra
# using ..Fundamentals # Předpokládá export BlockMesh, eval_sdf, quantize atd.

# Funkce cut_edge! zůstává stejná, protože identifikuje kladný/záporný vrchol
# a interpoluje mezi nimi bez ohledu na to, která hodnota znamená "uvnitř".
function cut_edge!(
    i::Int, j::Int,
    mesh::BlockMesh,
    node_sdf::Vector{Float64}, # Explicitní předání node_sdf pro jasnost
    cut_map::Dict{Tuple{Int, Int}, Int}
)::Int
    # ... (kód cut_edge! zůstává beze změny) ...
    # Zajištění, že SDF hodnoty mají opačná znaménka nebo je jedna nulová
    sdf_i = node_sdf[i]
    sdf_j = node_sdf[j]

    # Zpracování případů, kdy je jeden vrchol již na povrchu (nebo velmi blízko)
    tol = mesh.grid_tol # Použití malé tolerance pro kontrolu nuly
    if abs(sdf_i) < tol
        return i # Pokud je uzel i již na povrchu, "řez" je jen uzel i
    elseif abs(sdf_j) < tol
        return j # Pokud je uzel j již na povrchu, "řez" je jen uzel j
    end

    # Pokračovat pouze pokud jsou znaménka striktně opačná
    if sign(sdf_i) == sign(sdf_j)
        error("cut_edge! voláno s vrcholy na stejné straně izoplochy: i=$i (sdf=$(sdf_i)), j=$j (sdf=$(sdf_j))")
    end

    # Kanonická reprezentace hrany (seřazené indexy)
    edge = (min(i, j), max(i, j))

    # Kontrola, zda tato hrana již byla řezána
    if haskey(cut_map, edge)
        return cut_map[edge]
    end

    # Výpočet interpolačního faktoru
    pos_i = mesh.X[i]
    pos_j = mesh.X[j]

    # Nechť vrchol 'a' je ten s kladným SDF (uvnitř), 'b' ten se záporným SDF (vně)
    local pos_a::SVector{3, Float64}, pos_b::SVector{3, Float64}
    local sdf_a::Float64, sdf_b::Float64
    if sdf_i > 0 # i je kladný (uvnitř), j je záporný (vně)
        pos_a, sdf_a = pos_i, sdf_i
        pos_b, sdf_b = pos_j, sdf_j
    else # j je kladný (uvnitř), i je záporný (vně)
        pos_a, sdf_a = pos_j, sdf_j
        pos_b, sdf_b = pos_i, sdf_i
    end

    # Lineární interpolace pro nalezení průsečíku s nulou
    # t = sdf_a / (sdf_a - sdf_b) # Váha pro (pos_b - pos_a)
    # new_pos = pos_a + t * (pos_b - pos_a)
    alpha = sdf_a / (sdf_a - sdf_b) # C++ alpha odpovídá váze pos_j, pokud i je kladné
    if sdf_i > 0 # alpha je váha pro pos_j
         new_pos = (1.0 - alpha) * pos_i + alpha * pos_j
    else # alpha je váha pro pos_i (protože j je nyní kladné)
         new_pos = alpha * pos_i + (1.0 - alpha) * pos_j
    end


    # Kvantizace nové pozice pro kontrolu existujících uzlů velmi blízko
    p_key = quantize(new_pos, mesh.grid_tol)

    local new_index::Int
    if haskey(mesh.node_hash, p_key)
        # Znovu použít existující uzel
        new_index = mesh.node_hash[p_key]
        # Zajistit, že jeho SDF je označeno jako nula, pokud ho znovu používáme
        if abs(mesh.node_sdf[new_index]) > tol
             mesh.node_sdf[new_index] = 0.0
        end
    else
        # Vytvořit nový uzel
        push!(mesh.X, new_pos)
        push!(mesh.node_sdf, 0.0) # Nový vrchol je na izoploše
        new_index = length(mesh.X)
        mesh.node_hash[p_key] = new_index
        # Je třeba také změnit velikost INE, pokud přidáváme uzly během této fáze - předpokládáme, že INE bude přestavěno později
    end

    # Uložit výsledek do cut_map pro tuto hranu
    cut_map[edge] = new_index
    return new_index
end


# Funkce, která aplikuje logiku "trim spikes" na jeden tetraedr
# Vrací vektor nových tetraedrů nahrazujících vstupní.
# Modifikuje mesh a cut_map prostřednictvím cut_edge!
function apply_stencil_trim_spikes!(
    mesh::BlockMesh,
    tet::Vector{Int64},
    cut_map::Dict{Tuple{Int, Int}, Int}
)::Vector{Vector{Int64}}

    node_indices = tet
    node_sdf = [mesh.node_sdf[idx] for idx in node_indices]
    tol = mesh.grid_tol # Tolerance pro kontrolu nuly

    # --- Kontrola speciálních případů ---

    # 1. Případ: Všechny vrcholy jsou uvnitř nebo na hranici (SDF >= 0)
    if all(s -> s >= -tol, node_sdf)
        # Pokud jsou všechny vrcholy uvnitř nebo na hranici, tetraedr je zcela uvnitř.
        # Ponechat původní tetraedr.
        return [tet]
    end

    # 2. Případ: Všechny vrcholy jsou vně (SDF < 0)
    if all(s -> s < tol, node_sdf)
        # Pokud jsou všechny vrcholy striktně vně, tetraedr zahodit.
        return Vector{Vector{Int64}}()
    end

    # --- Obecný případ: Tetraedr protíná izoplochu ---

    # Vytvořit páry (SDF hodnota, původní index) pro třídění
    vert_data = [(node_sdf[i], node_indices[i]) for i in 1:4]

    # Třídit vrcholy: VZESTUPNĚ podle SDF, remízy řešit vzestupným indexem.
    # Nejnižší (nejvíce vně) SDF bude první.
    sort!(vert_data, by = x -> (x[1], x[2]))

    # Extrahovat seřazené indexy a jejich SDF hodnoty
    # s, r, q, p kde sdf_s <= sdf_r <= sdf_q <= sdf_p
    s, r, q, p = [vd[2] for vd in vert_data]
    sdf_s, sdf_r, sdf_q, sdf_p = [vd[1] for vd in vert_data]

    # Určit orientaci (flipped) - prozatím vynecháno, předpokládáme standardní
    flipped = false # Placeholder

    # --- Analýza případů podle znamének SDF (P=Positive/Inside, N=Negative/Outside, Z=Zero) ---
    is_s_neg = sdf_s < -tol
    is_r_neg = sdf_r < -tol
    is_q_neg = sdf_q < -tol
    is_p_neg = sdf_p < -tol # Mělo by být false, protože jsme již odfiltrovali případ NNNN

    is_s_zero = abs(sdf_s) <= tol
    is_r_zero = abs(sdf_r) <= tol
    is_q_zero = abs(sdf_q) <= tol
    is_p_zero = abs(sdf_p) <= tol

    # Kategorie vrcholů: P (kladný/uvnitř), N (záporný/vně), Z (nulový/na hranici)
    s_cat = is_s_neg ? 'N' : (is_s_zero ? 'Z' : 'P')
    r_cat = is_r_neg ? 'N' : (is_r_zero ? 'Z' : 'P')
    q_cat = is_q_neg ? 'N' : (is_q_zero ? 'Z' : 'P')
    p_cat = is_p_neg ? 'N' : (is_p_zero ? 'Z' : 'P') # Měl by být P nebo Z

    pattern = (s_cat, r_cat, q_cat, p_cat)
    # println("Pattern: ", pattern, " Indices: ", [s, r, q, p], " SDFs: ", [sdf_s, sdf_r, sdf_q, sdf_p])

    new_tets = Vector{Vector{Int64}}()

    # Cílem je vytvořit nové tetraedry, které pokrývají část původního tetraedru,
    # kde je SDF >= 0 (tj. část uvnitř nebo na hranici).

    if pattern == ('N', 'N', 'N', 'P') || pattern == ('N', 'N', 'N', 'Z') # ---+ (3 vně, 1 uvnitř/na hranici)
        # Vrchol p je uvnitř/na hranici.
        sp = cut_edge!(s, p, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r, p, mesh, mesh.node_sdf, cut_map)
        qp = cut_edge!(q, p, mesh, mesh.node_sdf, cut_map)
        # Nový tetraedr tvořený řeznými body a vnitřním vrcholem p
        push!(new_tets, [sp, rp, qp, p]) # TODO: Zkontrolovat orientaci ok
        # push!(new_tets, [sp, qp, rp, p]) # TODO: Zkontrolovat orientaci

    elseif pattern == ('N', 'N', 'P', 'P') || pattern == ('N', 'N', 'Z', 'P') || pattern == ('N', 'N', 'P', 'Z') || pattern == ('N', 'N', 'Z', 'Z') # --++ (2 vně, 2 uvnitř/na hranici)
        # Vrcholy q a p jsou uvnitř/na hranici.
        sq = cut_edge!(s, q, mesh, mesh.node_sdf, cut_map)
        sp = cut_edge!(s, p, mesh, mesh.node_sdf, cut_map)
        rq = cut_edge!(r, q, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r, p, mesh, mesh.node_sdf, cut_map)
        # Čtyřúhelník řezu: sq, sp, rp, rq. Základna: q, p.
        # Tetraedrizace hranolu (řez + základna). Rozdělení čtyřúhelníku úhlopříčkou
        # vycházející z řezu nejvíce "vnitřního" vrcholu (p). Úhlopříčka sp-rq.
        # Tet 1: [q, p, rq, sq] # Základna q, rq, sq
        # Tet 2: [p, rp, rq, sp] # Základna p, rp, rq - chyba v logice
        # Tet 3: [p, sq, sp, rq] # Základna p, sq, sp - chyba v logice

        # Použití C++ logiky ++-- s převrácenými rolemi P/N:
        # Původní: (pr, qr, r, s), (pr, qs, qr, s), (pr, ps, qs, s) (r,s byly uvnitř)
        # Nyní: p,q jsou uvnitř; s,r jsou vně. Řezy: sq, sp, rq, rp.
        # Chceme část s p,q.
        # push!(new_tets, [p, q, rq, rp]) # Tet 1 (základna p, q, rq)
        # push!(new_tets, [p, rq, sq, rp]) # Tet 2 (základna p, rq, sq)
        # push!(new_tets, [p, sq, sp, rp]) # Tet 3 (základna p, sq, sp)

        push!(new_tets, [rp, rq, q, p]) # ok?
        push!(new_tets, [rp, sq, rq, p])
        push!(new_tets, [rp, sp, sq, p])

        # push!(new_tets, [p, q, rq, sq])
        # push!(new_tets, [p, q, sp, sq])
        # push!(new_tets, [p, q, rp, sp])
        
        # push!(new_tets, [p, q, rp, sp])
        # push!(new_tets, [p, q, rq, rp])
        # push!(new_tets, [p, q, sq, rq])

        # TODO: Zkontrolovat orientaci a konzistenci triangulace

    elseif pattern == ('N', 'P', 'P', 'P') || pattern == ('N', 'Z', 'P', 'P') || pattern == ('N', 'P', 'Z', 'P') || pattern == ('N', 'P', 'P', 'Z') || pattern == ('N', 'Z', 'Z', 'P') || pattern == ('N', 'Z', 'P', 'Z') || pattern == ('N', 'P', 'Z', 'Z') || pattern == ('N', 'Z', 'Z', 'Z') # -+++ (1 vně, 3 uvnitř/na hranici)
        # Vrchol s je vně. Vrcholy r, q, p jsou uvnitř/na hranici.
        sr = cut_edge!(s, r, mesh, mesh.node_sdf, cut_map)
        sq = cut_edge!(s, q, mesh, mesh.node_sdf, cut_map)
        sp = cut_edge!(s, p, mesh, mesh.node_sdf, cut_map)
        # Trojúhelník řezu: sr, sq, sp. Vnitřní vrcholy: r, q, p.
        # Nahradit původní tetraedr třemi novými.
        # Použití C++ logiky +--- s převrácenými rolemi P/N:
        # Původní: (pq, q, r, s), (pq, r, pr, s), (pq, s, pr, ps) (p byl uvnitř)
        # Nyní: r,q,p jsou uvnitř; s je vně. Řezy: sr, sq, sp.
        # Chceme část s r,q,p.
        # push!(new_tets, [p, q, r, sr]) # Tet 1 (základna p, q, r)
        # push!(new_tets, [p, q, sr, sq]) # Tet 2 (základna p, q, sr)
        # push!(new_tets, [p, sr, sp, sq]) # Tet 3 (základna p, sr, sp)
        
        push!(new_tets, [sr, r, q, p]) # ok (pravděpodobně víc ok než NNPP)
        push!(new_tets, [sr, q, sq, p])
        push!(new_tets, [sr, p, sq, sp])

        # push!(new_tets, [p, r, q, sr])
        # push!(new_tets, [p, q, sq, sr])
        # push!(new_tets, [p, sp, sq, sr])
        
        # push!(new_tets, [sr, r, q, p])
        # push!(new_tets, [sq, sr, q, p])
        # push!(new_tets, [sp, sr, sq, p])

        # TODO: Zkontrolovat orientaci a konzistenci triangulace

    elseif pattern == ('Z', 'Z', 'Z', 'Z') # ZZZZ (všechny na hranici)
        # Speciální případ: tetraedr leží celý na izoploše.
        # Zkontrolovat centroid, zda je uvnitř nebo vně.
        centroid = (mesh.X[s] + mesh.X[r] + mesh.X[q] + mesh.X[p]) / 4.0
        if eval_sdf(mesh, centroid) < -tol # Centroid je vně
             return Vector{Vector{Int64}}() # Zahodit vnější povrchový tetraedr
        else
             # return [tet] # Ponechat vnitřní povrchový tetraedr ok
             return Vector{Vector{Int64}}() # Zahodit vnější povrchový tetraedr
        end

    else
        # Neočekávaný vzor - může nastat, pokud `cut_edge` selže nebo tolerance způsobí problémy.
        # Nebo pokud případ NNNN/PPPP nebyl správně odfiltrován na začátku.
        @warn "Neočekávaný vzor SDF v apply_stencil_trim_spikes!: $pattern. Tet: $tet. SDFs: $node_sdf"
        # Rozhodnutí: zahodit nebo ponechat? Bezpečnější je zahodit.
        return Vector{Vector{Int64}}()
    end

    # Volitelné: Zkontrolovat a opravit orientaci nových tetraedrů
    # for nt in new_tets
    #     check_and_fix_orientation!(nt, mesh)
    # end

    return new_tets
end

# (Zachovejte ostatní funkce jako process_cell_A15!, process_cell_Schlafli!, atd.)
