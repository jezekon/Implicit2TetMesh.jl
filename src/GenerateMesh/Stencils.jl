# Stencils.jl
# (Zachovejte existující includes a definice struktur)
# using StaticArrays
# using LinearAlgebra
# using ..Fundamentals # Předpokládá export BlockMesh, eval_sdf, quantize atd.

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



# Funkce cut_edge! zůstává stejná
function cut_edge!(
    i::Int, j::Int,
    mesh::BlockMesh,
    node_sdf::Vector{Float64}, # Explicitní předání node_sdf pro jasnost
    cut_map::Dict{Tuple{Int, Int}, Int}
)::Int
    # Zajištění, že SDF hodnoty mají opačná znaménka nebo je jedna nulová
    sdf_i = node_sdf[i]
    sdf_j = node_sdf[j]

    tol = mesh.grid_tol # Použití malé tolerance pro kontrolu nuly
    is_i_zero = abs(sdf_i) < tol
    is_j_zero = abs(sdf_j) < tol

    # Zpracování případů, kdy je jeden vrchol již na povrchu
    if is_i_zero && is_j_zero
        # Oba jsou na povrchu - technicky by se nemělo řezat, ale vrátíme jeden z nich.
        # Může indikovat degenerovaný případ nebo potřebu úpravy tolerance.
        # Vrácení vrcholu s menším indexem pro konzistenci.
        return min(i, j)
    elseif is_i_zero
        return i # Pokud je uzel i již na povrchu, "řez" je jen uzel i
    elseif is_j_zero
        return j # Pokud je uzel j již na povrchu, "řez" je jen uzel j
    end

    # Pokračovat pouze pokud jsou znaménka striktně opačná
    if sign(sdf_i) == sign(sdf_j)
        # Pokud mají stejné znaménko a ani jeden není nulový, je to chyba
        error("cut_edge! voláno s vrcholy na stejné straně izoplochy (a žádný není nulový): i=$i (sdf=$(sdf_i)), j=$j (sdf=$(sdf_j))")
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

    # Lineární interpolace pro nalezení průsečíku s nulou
    # alpha je váha pro pos_j, pokud sdf_i > 0
    # alpha je váha pro pos_i, pokud sdf_j > 0
    alpha = sdf_i / (sdf_i - sdf_j)
    # Zajistit, aby alpha byla v [0, 1] i při numerických nepřesnostech blízko 0
    alpha = clamp(alpha, 0.0, 1.0)

    # Interpolovaná pozice
    # new_pos = (1.0 - alpha) * pos_i + alpha * pos_j # Toto je správně obecně
    # Přesnější kontrola podle C++ logiky (i když výsledek by měl být stejný)
    if sdf_i > 0 # i je kladný (uvnitř), j je záporný (vně)
         # alpha je váha pro pos_j
         new_pos = (1.0 - alpha) * pos_i + alpha * pos_j
    else # j je kladný (uvnitř), i je záporný (vně)
         # alpha je váha pro pos_i
         # Vzorec sdf_i / (sdf_i - sdf_j) dává váhu pro druhý bod (j),
         # takže (1-alpha) je váha pro i.
         # Ale C++ kód počítá alpha = vphi[i]/(vphi[i]-vphi[j]) a pak
         # (1-alpha)*V(i) + alpha*V(j). Takže Julia kód je konzistentní.
         new_pos = (1.0 - alpha) * pos_i + alpha * pos_j
    end


    # Kvantizace nové pozice pro kontrolu existujících uzlů velmi blízko
    p_key = quantize(new_pos, mesh.grid_tol)

    local new_index::Int
    if haskey(mesh.node_hash, p_key)
        # Znovu použít existující uzel
        new_index = mesh.node_hash[p_key]
        # Zajistit, že jeho SDF je označeno jako nula, pokud ho znovu používáme
        # a není již nulové (kvůli toleranci)
        if abs(mesh.node_sdf[new_index]) > tol
             mesh.node_sdf[new_index] = 0.0
        end
    else
        # Vytvořit nový uzel
        push!(mesh.X, new_pos)
        push!(mesh.node_sdf, 0.0) # Nový vrchol je na izoploše
        new_index = length(mesh.X)
        mesh.node_hash[p_key] = new_index
        # Poznámka: INE bude třeba aktualizovat později hromadně
    end

    # Uložit výsledek do cut_map pro tuto hranu
    cut_map[edge] = new_index
    return new_index
end


# --- Zde vložte definice funkcí ---
# check_tetrahedron_orientation(mesh::BlockMesh, tet::Vector{Int})
# fix_tetrahedron_orientation!(mesh::BlockMesh, tet::Vector{Int})
# --- Konec definic funkcí ---

"""
    check_tetrahedron_orientation(mesh::BlockMesh, tet::Vector{Int})

Checks if a tetrahedron has positive orientation (positive Jacobian determinant).
Returns true for correctly oriented tetrahedra, false for inverted ones.
"""
function check_tetrahedron_orientation(mesh::BlockMesh, tet::Vector{Int})
    # Get vertices of the tetrahedron
    # Ensure we don't try to access index 0 or out of bounds
    if length(tet) != 4 || any(i -> i <= 0 || i > length(mesh.X), tet)
        @warn "Invalid tetrahedron indices for orientation check: $tet"
        return false # Consider invalid index as bad orientation
    end
    vertices = [mesh.X[tet[i]] for i in 1:4]

    # Calculate edge vectors from first vertex
    a = vertices[2] - vertices[1]
    b = vertices[3] - vertices[1]
    c = vertices[4] - vertices[1]

    # Calculate the Jacobian determinant (volume is det/6)
    # Use a small tolerance to avoid issues with near-zero volume tets
    det_value = dot(a, cross(b, c))

    # Return true for positive orientation, false for negative or zero
    # Consider zero volume as potentially problematic too, though not strictly inverted.
    # Adjust tolerance as needed. Using 1e-12 as a small positive threshold.
    return det_value > 1e-12
end

"""
    fix_tetrahedron_orientation!(mesh::BlockMesh, tet::Vector{Int})

Fixes the orientation of a tetrahedron by swapping two vertices if needed.
This modifies the input tetrahedron in-place. Returns true if fixed, false otherwise.
"""
function fix_tetrahedron_orientation!(mesh::BlockMesh, tet::Vector{Int})
    if length(tet) != 4
         @warn "Attempting to fix orientation of non-tetrahedron: $tet"
         return false
    end
    if !check_tetrahedron_orientation(mesh, tet)
        # Swap the last two vertices to change orientation
        # This is a common convention for flipping orientation.
        tet[3], tet[4] = tet[4], tet[3]
        # Optional: Double-check if the flip actually worked
        # if !check_tetrahedron_orientation(mesh, tet)
        #     @warn "Failed to fix orientation for tet: $tet after swapping 3 and 4."
        #     # Potentially try swapping other pairs? Or just report.
        # end
        return true  # Orientation was potentially fixed
    end
    return false # No fix needed
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
    # Získání SDF hodnot přímo z mesh.node_sdf pro konzistenci
    node_sdf = [mesh.node_sdf[idx] for idx in node_indices]
    tol = mesh.grid_tol # Tolerance pro kontrolu nuly

    # --- Kontrola speciálních případů ---
    if all(s -> s < -tol, node_sdf) # Všechny striktně záporné (vně)
        return Vector{Vector{Int64}}() # Zahodit
    end
    if all(s -> s >= -tol, node_sdf) # Všechny nezáporné (uvnitř nebo na hranici)
        # Ponechat původní tet, ale zkontrolovat jeho orientaci pro jistotu
        # (I když by původní tety měly být v pořádku)
        # fix_tetrahedron_orientation!(mesh, tet) # Volitelné zde
        return [tet] # Ponechat
    end
    # Nyní víme, že tet protíná hranici (některé < -tol, některé >= -tol)

    # --- Obecný případ: Tetraedr protíná izoplochu ---

    vert_data = [(node_sdf[i], node_indices[i]) for i in 1:4]
    p = [1, 2, 3, 4]
    flipped = false
    less_than(idx1, idx2) = (vert_data[idx1][1] < vert_data[idx2][1]) || (vert_data[idx1][1] == vert_data[idx2][1] && vert_data[idx1][2] < vert_data[idx2][2])

    if less_than(p[2], p[1]) p[1], p[2] = p[2], p[1]; flipped = !flipped end
    if less_than(p[4], p[3]) p[3], p[4] = p[4], p[3]; flipped = !flipped end
    if less_than(p[3], p[1]) p[1], p[3] = p[3], p[1]; flipped = !flipped end
    if less_than(p[4], p[2]) p[2], p[4] = p[4], p[2]; flipped = !flipped end
    if less_than(p[3], p[2]) p[2], p[3] = p[3], p[2]; flipped = !flipped end

    s_idx, r_idx, q_idx, p_idx = [vert_data[pi][2] for pi in p]
    sdf_s, sdf_r, sdf_q, sdf_p = [vert_data[pi][1] for pi in p]

    is_s_neg = sdf_s < -tol
    is_r_neg = sdf_r < -tol
    is_q_neg = sdf_q < -tol
    # is_p_neg = sdf_p < -tol # Mělo by být false

    is_p_pos = sdf_p > tol
    is_q_pos = sdf_q > tol
    is_r_pos = sdf_r > tol
    is_s_pos = sdf_s > tol # Mělo by být false

    is_s_zero = abs(sdf_s) <= tol
    is_r_zero = abs(sdf_r) <= tol
    is_q_zero = abs(sdf_q) <= tol
    is_p_zero = abs(sdf_p) <= tol

    new_tets = Vector{Vector{Int64}}()

    # --- Helper function to add tet with orientation check ---
    function add_tet!(t::Vector{Int})
        # Check for duplicate nodes within the tet before fixing orientation
        if length(Set(t)) != 4
             @warn "Degenerate tetrahedron generated (duplicate nodes): $t. Skipping."
             return
        end
        fix_tetrahedron_orientation!(mesh, t)
        # Optional: Add a final check after fixing
        if !check_tetrahedron_orientation(mesh, t)
             @warn "Tetrahedron $t still has incorrect/zero orientation after fix attempt. Skipping."
             return
        end
        push!(new_tets, t)
    end

    # --- Implementace případů podle C++ logiky s opravou orientace ---

    if is_s_zero && !is_r_neg && !is_q_neg && !is_p_neg # Ekvivalent C++ +++0
        # Zpracování ZZZZ a ostatních případů jako předtím
        if is_r_zero && is_q_zero && is_p_zero # ZZZZ
             centroid = (mesh.X[s_idx] + mesh.X[r_idx] + mesh.X[q_idx] + mesh.X[p_idx]) / 4.0
             if eval_sdf(mesh, centroid) < -tol # Centroid je vně
                  return Vector{Vector{Int64}}()
             else
                  return Vector{Vector{Int64}}() # Zahodit povrchové tety
             end
        else # ZPPP, ZZPP, ZZZP
             return Vector{Vector{Int64}}() # Zahodit
        end

    elseif is_s_neg && is_r_neg && is_q_neg && (is_p_pos || is_p_zero) # NNNP / NNNZ (Ekvivalent C++ +++-)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        qp = cut_edge!(q_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        if flipped
            add_tet!([rp, sp, qp, p_idx])
        else
            add_tet!([sp, rp, qp, p_idx])
        end

    elseif is_s_neg && (is_p_pos || is_p_zero) && (is_q_pos || is_q_zero) && (is_r_pos || is_r_zero) # NPPP / NPZZ / atd. (Ekvivalent C++ +---)
        # Tento případ by měl být is_s_neg a ostatní >= -tol
        # Kontrola: !(is_r_neg || is_q_neg || is_p_neg)
        if !is_r_neg && !is_q_neg # Potvrzení, že r, q, p jsou >= -tol
            sr = cut_edge!(s_idx, r_idx, mesh, mesh.node_sdf, cut_map)
            sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
            sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
            # Tři tetraedry tvořící vnitřní část
            if flipped
                add_tet!([q_idx, r_idx, p_idx, sr])
                add_tet!([p_idx, q_idx, sr, sq])
                add_tet!([sr, p_idx, sq, sp])
            else
                add_tet!([r_idx, q_idx, p_idx, sr])
                add_tet!([q_idx, p_idx, sr, sq])
                add_tet!([p_idx, sr, sq, sp])
            end
        else
             @warn "Logická chyba v NPPP větvi: r=$sdf_r, q=$sdf_q, p=$sdf_p"
             return Vector{Vector{Int64}}()
        end


    elseif is_s_neg && is_r_neg && (is_q_pos || is_q_zero) && (is_p_pos || is_p_zero) # NNPP / NNZP / NNPZ / NNZZ (Ekvivalent C++ ++--)
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rq = cut_edge!(r_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        # Šest tetraedrů tvořících vnitřní část (triangulace čtyřstěnu)
        # Použijeme C++ triangulaci (přeloženou) a opravíme orientaci
        # C++ Tety: (pr, qr, r, s), (pr, qs, qr, s), (pr, ps, qs, s)
        # Julia překlad (s,r záporné; q,p kladné):
        # Tet 1: (sq, rq, q, p) ? Ne, to jsou 4 kladné/nulové body.
        # Tet 1: (q, p, rq, rp)
        # Tet 2: (q, p, rp, sp) ? Ne.
        # Tet 2: (q, p, sp, sq) ? Ne.
        # Tet 3: (q, rq, rp, sp) ? Ne.
        # Tet 3: (q, rq, sp, sq) ? Ne.
        # Zkusme triangulaci s vrcholy q, p a řezy sq, sp, rq, rp
        # Tvoříme tety s vrcholy q a p.
        # Tet 1: (q, p, rq, rp) # Společná stěna q-p-rq
        # Tet 2: (q, p, rp, sp) # Společná stěna q-p-rp
        # Tet 3: (q, p, sp, sq) # Společná stěna q-p-sp
        # Toto dává 3 tety. C++ má 3 tety. Zkusme to.
        # C++ (pr, qr, r, s) -> (sq, rq, q_idx, p_idx) ? Ne.
        # C++ (pr, qs, qr, s) -> (sq, sp, rq, p_idx) ? Ne.
        # C++ (pr, ps, qs, s) -> (sq, rp, sp, p_idx) ? Ne.

        # Použijeme C++ triangulaci z kódu (i když se zdá divná) a opravíme
        if flipped
            add_tet!([p_idx, rq, q_idx, sq]) # Tet 1 flipped
            add_tet!([p_idx, sp, sq, rp])    # Tet 2 flipped
            add_tet!([p_idx, rp, sq, rq])    # Tet 3 flipped
        else
            add_tet!([p_idx, q_idx, rq, sq]) # Tet 1
            add_tet!([p_idx, sq, sp, rp])    # Tet 2
            add_tet!([p_idx, sq, rp, rq])    # Tet 3
        end

    # --- Zpracování případů s nulami ---

    elseif is_s_neg && is_r_zero && (is_q_pos || is_q_zero) && (is_p_pos || is_p_zero) # NZPP / NZPZ / NZZP / NZZZ (Ekvivalent C++ ++0-)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        # C++ Tet: (ps, qs, r, s) -> Julia (sp, sq, r_idx, p_idx)? Ne.
        # Logika: Vnitřní část má vrcholy r, q, p a řezy sp, sq.
        # Tvoříme tety s vrcholy r, q, p.
        # Tet 1: (r_idx, q_idx, p_idx, sq)
        # Tet 2: (r_idx, p_idx, sq, sp)
        if flipped
             add_tet!([q_idx, r_idx, p_idx, sq]) # Tet 1 flipped?
             add_tet!([p_idx, r_idx, sq, sp]) # Tet 2 flipped?
        else
             add_tet!([r_idx, q_idx, p_idx, sq]) # Tet 1
             add_tet!([r_idx, p_idx, sq, sp]) # Tet 2
        end


    elseif is_s_neg && is_r_neg && is_q_zero && (is_p_pos || is_p_zero) # NNZP / NNZZ (Ekvivalent C++ +0--)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        # C++ Tets: (pr, q, r, s), (ps, q, pr, s)
        # Julia překlad: (rp, q_idx, r_idx, s_idx), (sp, q_idx, rp, s_idx) ? Ne.
        # Logika: Vnitřní část má vrcholy q, p a řezy sp, rp.
        # Tet 1: (q_idx, p_idx, rp, sp)
        if flipped
             add_tet!([p_idx, q_idx, rp, sp]) # Tet 1 flipped?
        else
             add_tet!([q_idx, p_idx, rp, sp]) # Tet 1
        end


    elseif is_s_neg && is_r_zero && is_q_zero && (is_p_pos || is_p_zero) # NZZP / NZZZ (Ekvivalent C++ +00-)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        # C++ Tet: (ps, q, r, s) -> Julia (sp, q_idx, r_idx, s_idx)? Ne.
        # Logika: Vnitřní část má vrcholy r, q, p a řez sp.
        # Tet 1: (r_idx, q_idx, p_idx, sp)
        if flipped
             add_tet!([q_idx, r_idx, p_idx, sp]) # Tet 1 flipped?
        else
             add_tet!([r_idx, q_idx, p_idx, sp]) # Tet 1
        end

    else
        @warn "Neočekávaný vzor SDF v apply_stencil_trim_spikes! Po třídění: s=$sdf_s, r=$sdf_r, q=$sdf_q, p=$sdf_p. Indices: s=$s_idx, r=$r_idx, q=$q_idx, p=$p_idx. Tet: $tet. Flipped: $flipped"
        return Vector{Vector{Int64}}() # Zahodit neočekávaný případ
    end

    return new_tets
end

# --- Zbytek kódu (slice_ambiguous_tetrahedra!, cut_edge!, atd.) zůstává ---
# Nezapomeňte zahrnout definice check_tetrahedron_orientation a fix_tetrahedron_orientation!
# a také definici BlockMesh a souvisejících struktur/funkcí.
