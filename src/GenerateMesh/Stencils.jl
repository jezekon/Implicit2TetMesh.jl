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

    # 1. Případ: Všechny vrcholy jsou uvnitř nebo na hranici (SDF >= -tol)
    #    V C++ je to vphi[p]<=0 && vphi[q]<=0 && vphi[r]<=0 && vphi[s]<=0
    #    což znamená, že tetraedr *nevyčnívá* (je celý uvnitř nebo na hranici)
    #    V Julia konvenci (SDF > 0 uvnitř) to znamená, že všechny SDF jsou <= tol
    #    ALE C++ kód zde kontroluje, zda *nevyčnívá*, tj. zda není potřeba řezat.
    #    Řezat je potřeba, pokud je *alespoň jeden* vrchol kladný (vně v C++, uvnitř v Julii).
    #    Takže C++ `continue` znamená "ponechat tet, pokud jsou všechny <= 0".
    #    Julia ekvivalent: Ponechat tet, pokud jsou všechny SDF <= tol.
    # if all(s -> s <= tol, node_sdf) # Původní C++ logika: pokud nic nevyčnívá (vše <= 0), nedělej nic
    #     # Tento tet je celý uvnitř nebo na hranici, není třeba řezat.
    #     # Ale my chceme jen část S >= 0. Takže pokud je celý S <= 0, měl by být zahozen?
    #     # Ne, C++ kód ho ponechává. Zdá se, že cílem C++ je odstranit části > 0.
    #     # Cílem Julia kódu je ponechat části >= 0.
    #     # Takže pokud jsou všechny <= tol, ale ne všechny < -tol, měl by být ponechán?
    #     # Pokud jsou všechny < -tol, zahodit.
    #     # Pokud jsou některé mezi [-tol, tol] a zbytek < -tol, ponechat.
    #     # Pokud jsou všechny <= tol, ponechat.
    #     # Zkusme se držet Julia cíle: ponechat část >= 0.
    #     if all(s -> s < -tol, node_sdf) # Všechny striktně záporné (vně)
    #          return Vector{Vector{Int64}}() # Zahodit
    #     else # Alespoň jeden je >= -tol (na hranici nebo uvnitř) a žádný není > tol
    #          return [tet] # Ponechat
    #     end
    # end
    # Jednodušší kontrola na začátku:
    if all(s -> s < -tol, node_sdf) # Všechny striktně záporné (vně)
        return Vector{Vector{Int64}}() # Zahodit
    end
    if all(s -> s >= -tol, node_sdf) # Všechny nezáporné (uvnitř nebo na hranici)
        return [tet] # Ponechat
    end
    # Nyní víme, že tet protíná hranici (některé < -tol, některé >= -tol)

    # --- Obecný případ: Tetraedr protíná izoplochu ---

    # Vytvořit páry (SDF hodnota, původní index) pro třídění
    vert_data = [(node_sdf[i], node_indices[i]) for i in 1:4]

    # Třídit vrcholy: VZESTUPNĚ podle SDF, remízy řešit vzestupným indexem.
    # Nejnižší (nejvíce vně) SDF bude první.
    # Implementace třídící sítě pro sledování `flipped` (jako v C++)
    # p0, p1, p2, p3 jsou indexy do vert_data (0-based pro snazší přepis)
    p = [1, 2, 3, 4] # Indexy do vert_data (1-based v Julii)
    flipped = false

    # Funkce pro porovnání (SDF, index)
    less_than(idx1, idx2) = (vert_data[idx1][1] < vert_data[idx2][1]) || (vert_data[idx1][1] == vert_data[idx2][1] && vert_data[idx1][2] < vert_data[idx2][2])

    # Třídící síť (přizpůsobená pro vzestupné třídění)
    if less_than(p[2], p[1]) p[1], p[2] = p[2], p[1]; flipped = !flipped end
    if less_than(p[4], p[3]) p[3], p[4] = p[4], p[3]; flipped = !flipped end
    if less_than(p[3], p[1]) p[1], p[3] = p[3], p[1]; flipped = !flipped end
    if less_than(p[4], p[2]) p[2], p[4] = p[4], p[2]; flipped = !flipped end
    if less_than(p[3], p[2]) p[2], p[3] = p[3], p[2]; flipped = !flipped end

    # Extrahovat seřazené indexy a jejich SDF hodnoty
    # s, r, q, p kde sdf_s <= sdf_r <= sdf_q <= sdf_p
    s_idx, r_idx, q_idx, p_idx = [vert_data[pi][2] for pi in p]
    sdf_s, sdf_r, sdf_q, sdf_p = [vert_data[pi][1] for pi in p]

    # --- Analýza případů podle znamének SDF ---
    # Použijeme přímé porovnání s tolerancí, abychom přesněji odpovídali C++
    is_s_neg = sdf_s < -tol
    is_r_neg = sdf_r < -tol
    is_q_neg = sdf_q < -tol
    is_p_neg = sdf_p < -tol # Mělo by být false, protože případ NNNN byl odfiltrován

    is_s_pos = sdf_s > tol
    is_r_pos = sdf_r > tol
    is_q_pos = sdf_q > tol
    is_p_pos = sdf_p > tol # Měl by být true, protože případ NNNN byl odfiltrován

    is_s_zero = abs(sdf_s) <= tol
    is_r_zero = abs(sdf_r) <= tol
    is_q_zero = abs(sdf_q) <= tol
    is_p_zero = abs(sdf_p) <= tol

    new_tets = Vector{Vector{Int64}}()

    # Přemapování na C++ logiku (kde p,q,r,s jsou seřazeny sestupně, p je nejkladnější)
    # Julia: s <= r <= q <= p (vzestupně)
    # C++:   p >= q >= r >= s (sestupně)
    # C++ `vphi[s] < 0` => Julia `sdf_s < -tol` (NNNN - zahozeno)
    # C++ `vphi[s] == 0` => Julia `is_s_zero` (a ostatní >= 0) => `+++0` (C++) => Zahodit tet
    # C++ `vphi[r] > 0` => Julia `sdf_r > tol` (a `sdf_s < -tol`) => `+++-` (C++) => Julia `NNNP`
    # C++ `vphi[q] < 0` => Julia `sdf_q < -tol` (a `sdf_p > tol`) => `+---` (C++) => Julia `NPPP`
    # C++ `vphi[q] > 0 && vphi[r] < 0` => Julia `sdf_q > tol` a `sdf_r < -tol` => `++--` (C++) => Julia `NNPP`
    # C++ `vphi[q] == 0 && vphi[r] == 0` => Julia `is_q_zero` a `is_r_zero` (a `sdf_s < -tol`, `sdf_p > tol`) => `+00-` (C++) => Julia `NZZP`
    # C++ `vphi[q] == 0` (a `vphi[r] < -tol`) => Julia `is_q_zero` a `sdf_r < -tol` (a `sdf_p > tol`) => `+0--` (C++) => Julia `NNZP`
    # C++ `vphi[r] == 0` (a `vphi[q] > tol`) => Julia `is_r_zero` a `sdf_q > tol` (a `sdf_s < -tol`) => `++0-` (C++) => Julia `NZZP` ? Ne, `NZPP`.

    # --- Implementace případů podle C++ logiky ---

    if is_s_zero && !is_r_neg && !is_q_neg && !is_p_neg # Ekvivalent C++ +++0 (všechny >= 0, nejnižší je 0)
        # Zkontrolovat centroid pro případ ZZZZ
        if is_r_zero && is_q_zero && is_p_zero # ZZZZ
             centroid = (mesh.X[s_idx] + mesh.X[r_idx] + mesh.X[q_idx] + mesh.X[p_idx]) / 4.0
             if eval_sdf(mesh, centroid) < -tol # Centroid je vně
                  return Vector{Vector{Int64}}() # Zahodit vnější povrchový tetraedr
             else
                  # Měli bychom ho ponechat? C++ `remove_exterior_tets` ho zahazuje, pokud je centroid > 0.
                  # Zdá se, že cílem je zahodit *všechny* tety ležící PŘESNĚ na povrchu.
                  return Vector{Vector{Int64}}() # Zahodit
             end
        else # Případy ZPPP, ZZPP, ZZZP - odpovídá C++ +++0
            # Zahodit tetraedr, který leží na hranici z vnější strany
             return Vector{Vector{Int64}}()
        end

    elseif is_s_neg && is_r_neg && is_q_neg && (is_p_pos || is_p_zero) # NNNP / NNNZ (Ekvivalent C++ +++-)
        # p_idx je uvnitř/na hranici, s_idx, r_idx, q_idx jsou vně
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        qp = cut_edge!(q_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        # C++: (ps, qs, rs, s) nebo (qs, ps, rs, s) [flipped]
        # Julia: p_idx odpovídá C++ s (záporný), {s,r,q}_idx odpovídá C++ {p,q,r} (kladné)
        #        sp, rp, qp odpovídá C++ ps, qs, rs
        # Nový tet: [sp, rp, qp, p_idx]
        if flipped
            push!(new_tets, [rp, sp, qp, p_idx]) # C++ (qs, ps, rs, s) -> Julia (rp, sp, qp, p_idx)
        else
            push!(new_tets, [sp, rp, qp, p_idx]) # C++ (ps, qs, rs, s) -> Julia (sp, rp, qp, p_idx)
        end

    elseif is_s_neg && (is_r_pos || is_r_zero) # NPPP / NPZZ / NPZP / NPPZ (Ekvivalent C++ +---)
        # s_idx je vně, r_idx, q_idx, p_idx jsou uvnitř/na hranici
        sr = cut_edge!(s_idx, r_idx, mesh, mesh.node_sdf, cut_map)
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        # C++: (pq, q, r, s), (pq, r, pr, s), (pq, s, pr, ps) nebo flipped
        # Julia: s_idx odpovídá C++ p (kladný), {r,q,p}_idx odpovídá C++ {q,r,s} (záporné)
        #        sr, sq, sp odpovídá C++ pq, pr, ps
        # Tet 1: [r_idx, q_idx, p_idx, sr] (základna r,q,p)
        # Tet 2: [q_idx, p_idx, sr, sq] (základna q,p,sr)
        # Tet 3: [p_idx, sr, sq, sp] (základna p,sr,sq)
        # Překlad C++ tetů do Julia indexů:
        # C++ (pq, q, r, s) -> Julia (sr, r_idx, q_idx, p_idx) ? Ne, to není správně.
        # C++ p je kladný vrchol, q,r,s záporné. pq, pr, ps jsou řezy.
        # C++ tety: (pq,q,r,s), (pq,r,pr,s), (pq,s,pr,ps)
        # Julia s je záporný, r,q,p kladné. sr, sq, sp jsou řezy.
        # Chceme část obsahující r,q,p. Trojúhelník řezu je sr,sq,sp.
        # Tet 1: (r, q, p, sr) - základna r,q,p; vrchol sr
        # Tet 2: (q, p, sr, sq) - základna q,p,sr; vrchol sq
        # Tet 3: (p, sr, sq, sp) - základna p,sr,sq; vrchol sp
        if flipped
            # Nutno prohodit druhé dva indexy v každém tetu pro zachování orientace
            push!(new_tets, [q_idx, r_idx, p_idx, sr]) # Tet 1 flipped
            push!(new_tets, [p_idx, q_idx, sr, sq]) # Tet 2 flipped
            push!(new_tets, [sr, p_idx, sq, sp]) # Tet 3 flipped
        else
            push!(new_tets, [r_idx, q_idx, p_idx, sr]) # Tet 1
            push!(new_tets, [q_idx, p_idx, sr, sq]) # Tet 2
            push!(new_tets, [p_idx, sr, sq, sp]) # Tet 3
        end


    elseif is_s_neg && is_r_neg && (is_q_pos || is_q_zero) # NNPP / NNZP / NNPZ / NNZZ (Ekvivalent C++ ++--)
        # s_idx, r_idx jsou vně, q_idx, p_idx jsou uvnitř/na hranici
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rq = cut_edge!(r_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        # C++: (pr, qr, r, s), (pr, qs, qr, s), (pr, ps, qs, s) nebo flipped
        # Julia: s,r odpovídá C++ p,q (kladné), q,p odpovídá C++ r,s (záporné)
        #        sq, sp, rq, rp odpovídá C++ pr, ps, qr, qs
        # Chceme část obsahující q_idx, p_idx. Čtyřúhelník řezu: sq, sp, rp, rq.
        # Triangulace C++:
        # Tet 1: (pr, qr, r, s) -> (sq, rq, q_idx, p_idx) ? Ne.
        # C++ p,q kladné; r,s záporné. Řezy pr, ps, qr, qs.
        # Tet 1: (pr, qr, r, s)
        # Tet 2: (pr, qs, qr, s)
        # Tet 3: (pr, ps, qs, s)
        # Julia s,r záporné; q,p kladné. Řezy sq, sp, rq, rp.
        # Chceme část s q, p.
        # Tet 1: (q, p, rq, rp) - základna q,p,rq; vrchol rp
        # Tet 2: (q, rp, sq, sp) - základna q,rp,sq; vrchol sp ? Ne.
        # Tet 2: (q, p, rp, sp) - základna q,p,rp; vrchol sp ? Ne.
        # Tet 2: (q, p, sp, sq) - základna q,p,sp; vrchol sq ? Ne.
        # Zkusme triangulovat čtyřúhelník řezu (sq, sp, rp, rq) úhlopříčkou sp-rq.
        # Tet 1: (q, p, rq, sp) # Základna q,p,rq
        # Tet 2: (q, sp, rq, sq) # Základna q,sp,rq
        # Tet 3: (p, sp, rp, rq) # Základna p,sp,rp
        # Porovnání s C++ tety (překlad indexů):
        # C++ Tet 1: (pr, qr, r, s) -> Julia (sq, rq, q, p)
        # C++ Tet 2: (pr, qs, qr, s) -> Julia (sq, sp, rq, p)
        # C++ Tet 3: (pr, ps, qs, s) -> Julia (sq, rp, sp, p)
        # Zdá se, že C++ triangulace je jiná. Použijme ji.
    # TODO: upravit popis podle Implementace níže
        if flipped
            # Prohodit druhé dva indexy
            push!(new_tets, [p_idx, rq, q_idx, sq]) # Tet 1 flipped
            push!(new_tets, [p_idx, sp, sq, rp])    # Tet 2 flipped
            push!(new_tets, [p_idx, rp, sq, rq])    # Tet 3 flipped
        else
            # Původní pořadí
            push!(new_tets, [p_idx, q_idx, rq, sq]) # Tet 1
            push!(new_tets, [p_idx, sq, sp, rp])    # Tet 2
            push!(new_tets, [p_idx, sq, rp, rq])    # Tet 3
        end
    # --- Zpracování případů s nulami, které nebyly pokryty výše ---
    # Tyto případy by měly odpovídat C++ ++0-, +00-, +0--

    elseif is_s_neg && is_r_zero && (is_q_pos || is_q_zero) # NZPP / NZPZ / NZZP / NZZZ (Ekvivalent C++ ++0-)
        # s_idx je vně, r_idx je na hranici, q_idx, p_idx jsou uvnitř/na hranici
        # C++: ++0- (p,q > 0; r = 0; s < 0). Řezy ps, qs. Tet (ps, qs, r, s)
        # Julia: s < 0; r = 0; q,p >= 0. Řezy sp, sq. Tet (sp, sq, r_idx, p_idx) ?
        # Překlad C++: ps -> sp, qs -> sq, r -> r_idx, s -> p_idx ? Ne.
        # C++ p,q kladné; r=0; s záporné.
        # Julia s záporné; r=0; q,p kladné.
        # C++ (ps, qs, r, s) -> Julia (sp, sq, r_idx, p_idx) ? Ne.
        # C++ s je záporný vrchol -> Julia p_idx
        # C++ r je nulový vrchol -> Julia r_idx
        # C++ p,q jsou kladné -> Julia q_idx, s_idx ? Ne, s_idx je záporný. q_idx, p_idx.
        # C++ ps, qs -> Julia řezy s kladnými: sq, sp.
        # C++ (ps, qs, r, s) -> Julia (sp, sq, r_idx, p_idx) ? Stále divné.
        # Zkusme logicky: Chceme část s r,q,p. Řezy jsou sp, sq. Základna r,q,p.
        # Tet (sp, sq, r_idx, p_idx) ? Ne.
        # Tet (r_idx, q_idx, p_idx, sp) ? Ne.
        # Tet (r_idx, q_idx, p_idx, sq) ? Ne.
        # Tet (q_idx, p_idx, r_idx, sq) ?
        # Tet (q_idx, p_idx, r_idx, sp) ?
        # C++ tet (ps, qs, r, s) má základnu (ps, qs, r) a vrchol s.
        # Julia ekvivalent: základna (sp, sq, r_idx), vrchol p_idx?
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        sq = cut_edge!(s_idx, q_idx, mesh, mesh.node_sdf, cut_map)
        if flipped
            push!(new_tets, [sq, sp, r_idx, p_idx]) # C++ (qs, ps, r, s)
        else
            push!(new_tets, [sp, sq, r_idx, p_idx]) # C++ (ps, qs, r, s)
        end


    elseif is_s_neg && is_r_neg && is_q_zero && (is_p_pos || is_p_zero) # NNZP / NNZZ (Ekvivalent C++ +0--)
        # s,r jsou vně, q je na hranici, p je uvnitř/na hranici
        # C++: +0-- (p > 0; q = 0; r,s < 0). Řezy pr, ps. Tets (pr, q, r, s), (ps, q, pr, s)
        # Julia: s,r < 0; q = 0; p >= 0. Řezy sp, rp.
        # C++ (pr, q, r, s) -> Julia (rp, q_idx, r_idx, p_idx) ? Ne.
        # C++ (ps, q, pr, s) -> Julia (sp, q_idx, rp, p_idx) ? Ne.
        # Logika: Chceme část s q,p. Řezy sp, rp.
        # Tet 1: (q_idx, p_idx, rp, r_idx) ?
        # Tet 2: (q_idx, p_idx, sp, rp) ?
        # Překlad C++:
        # Tet 1: (pr, q, r, s) -> (rp, q_idx, r_idx, s_idx) ?
        # Tet 2: (ps, q, pr, s) -> (sp, q_idx, rp, s_idx) ?
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        rp = cut_edge!(r_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        if flipped
            # Prohodit druhé dva indexy
            push!(new_tets, [q_idx, rp, r_idx, s_idx]) # Tet 1 flipped
            push!(new_tets, [q_idx, sp, rp, s_idx]) # Tet 2 flipped
        else
            push!(new_tets, [rp, q_idx, r_idx, s_idx]) # Tet 1
            push!(new_tets, [sp, q_idx, rp, s_idx]) # Tet 2
        end

    elseif is_s_neg && is_r_zero && is_q_zero && (is_p_pos || is_p_zero) # NZZP / NZZZ (Ekvivalent C++ +00-)
        # s je vně, r,q jsou na hranici, p je uvnitř/na hranici
        # C++: +00- (p > 0; q=0, r=0; s < 0). Řez ps. Tet (ps, q, r, s)
        # Julia: s < 0; r=0, q=0; p >= 0. Řez sp.
        # C++ (ps, q, r, s) -> Julia (sp, q_idx, r_idx, s_idx) ?
        sp = cut_edge!(s_idx, p_idx, mesh, mesh.node_sdf, cut_map)
        if flipped
            push!(new_tets, [q_idx, sp, r_idx, s_idx]) # C++ (ps, q, r, s) -> prohozeno q, sp
        else
            push!(new_tets, [sp, q_idx, r_idx, s_idx]) # C++ (ps, q, r, s)
        end

    else
        # Neočekávaný vzor - může nastat, pokud `cut_edge` selže nebo tolerance způsobí problémy.
        # Nebo pokud případ NNNN/PPPP nebyl správně odfiltrován na začátku.
        @warn "Neočekávaný vzor SDF v apply_stencil_trim_spikes! Po třídění: s=$sdf_s, r=$sdf_r, q=$sdf_q, p=$sdf_p. Indices: s=$s_idx, r=$r_idx, q=$q_idx, p=$p_idx. Tet: $tet. Flipped: $flipped"
        # Rozhodnutí: zahodit nebo ponechat? Bezpečnější je zahodit.
        return Vector{Vector{Int64}}()
    end

    # Volitelné: Zkontrolovat a opravit orientaci nových tetraedrů
    # V C++ kódu se orientace řeší pomocí `flipped`. Měli bychom to zde také udělat.
    # Kontrola orientace by měla být provedena po všech řezech.

    return new_tets
end

# (Zachovejte ostatní funkce jako process_cell_A15!, process_cell_Schlafli!, atd.)
