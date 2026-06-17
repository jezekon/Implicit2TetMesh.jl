# ==============================================================================
# Cap recovery ("warp #2") -- optional PRE-SLICE pass that reclaims inward spikes.
# DEFAULT OFF.
# ==============================================================================
# After warp! a thin layer of "+000 spike" tetrahedra survives: three vertices were
# warped onto phi = 0, the fourth (the apex) sits just OUTSIDE, yet the tet's centroid is
# INSIDE the material. The trimmer discards these (no vertex strictly inside, one outside),
# so the meshed solid loses a sliver of volume exactly where the surface pokes inward.
#
# This pass closes that gap WITHOUT touching connectivity. For each such spike it moves the
# shared outer apex onto phi = 0 along the average inward base normal and relabels it
# on-surface (node_sdf = 0) -- the exact same primitive warp! uses (move a node onto the
# surface, set its node_sdf to zero). The spike then becomes a "quadruple-zero" tet that the
# existing slice keeps via resolve_surface_candidates! (its centroid, now even further
# inside, tests inside; Labelle's [min,max] dihedral filter still rejects the few too-flat
# ones). Because it only ever sets node_sdf to 0 and moves the node onto phi = 0 -- the same
# state warp! produces -- the slice sees nothing it does not already handle, so the
# watertight cut structure (consistent vert sort + per-edge cut_map) is preserved by the
# same construction as the warp.
#
# SCOPE (v1): only the pure +000 pattern (exactly one outside vertex). The ambiguous ++00 /
# +++0 patterns (which of the several outside vertices is "the apex"?) are left for later;
# fans whose base normals diverge too much are skipped (the average normal is meaningless
# there) and accepted as residual.
#
# DETERMINISTIC: every apex is a DISTINCT node and is never another spike's base (an apex
# has node_sdf > tol, a base has node_sdf = 0, so the sets cannot overlap). The moves are
# therefore mutually independent -- they are all computed from the pre-move positions and
# applied at once, so the result does not depend on iteration order.

"""
    recover_caps!(mesh::BlockMesh; normal_spread_max_deg = 75.0) -> Int

Pre-slice cap-recovery pass (default OFF, see [`MeshGenerationOptions`](@ref)). Detects pure
"+000 spike" tetrahedra -- three vertices on `phi = 0`, exactly one outside, centroid inside
-- and snaps each spike's shared outer apex onto `phi = 0` along the average inward base
normal, relabelling it on-surface. The spike is thereby converted into a quadruple-zero
tetrahedron that the slice retains, recovering the inward volume the trimmer would discard.

Only `mesh.X` and `mesh.node_sdf` of the apex nodes change; connectivity is untouched, so it
must run AFTER warp's `update_connectivity!` and BEFORE `slice_ambiguous_tetrahedra!`
(the slice rebuilds its own cut structure from `node_sdf`, so no connectivity refresh is
needed in between). One O(N) pass over `mesh.IEN`; no inverse connectivity (`mesh.INE`)
required.

# Arguments
- `mesh::BlockMesh`: the warped lattice (post-`warp!`, pre-slice).

# Keywords
- `normal_spread_max_deg::Float64 = 75.0`: an apex is skipped (accepted as residual) when
  any two of its base-face normals diverge by more than this angle -- their average is not a
  meaningful move direction for such a spread-out fan.

# Returns
- the number of apexes actually moved onto `phi = 0`.
"""
function recover_caps!(mesh::BlockMesh; normal_spread_max_deg::Float64 = 75.0)
    isempty(mesh.IEN) && return 0
    tol = mesh.grid_tol   # SAME tolerance the slice classifies with, so the +000 set matches

    # (1)+(2) Detect pure +000 spikes with an inside centroid and group their on-surface base
    # triples by the shared outside apex. The classification mirrors the slice exactly:
    #   outside : node_sdf >  tol ,  inside : node_sdf < -tol ,  surface : |node_sdf| <= tol.
    bases = Dict{Int,Vector{NTuple{3,Int}}}()
    for tet in mesh.IEN
        no = ni = nz = 0
        apex = 0
        @inbounds for v in tet
            s = mesh.node_sdf[v]
            if s > tol
                no += 1
                apex = v
            elseif s < -tol
                ni += 1
            else
                nz += 1
            end
        end
        (no == 1 && ni == 0 && nz == 3) || continue

        # Keep only spikes that genuinely poke INTO the material (quartet's centroid test,
        # the same detector validated in scratch/volume_gap.jl).
        centroid =
            (mesh.X[tet[1]] + mesh.X[tet[2]] + mesh.X[tet[3]] + mesh.X[tet[4]]) / 4.0
        eval_sdf(mesh, centroid) < 0.0 || continue

        # The three base nodes are the tet's vertices other than the apex (tet order).
        b1 = b2 = b3 = 0
        @inbounds for v in tet
            v == apex && continue
            if b1 == 0
                b1 = v
            elseif b2 == 0
                b2 = v
            else
                b3 = v
            end
        end
        push!(get!(bases, apex, NTuple{3,Int}[]), (b1, b2, b3))
    end
    isempty(bases) && return 0

    # (3) For each apex: average inward base normal, fan-spread guard, re-project onto phi = 0.
    # Stage the accepted moves; do not apply yet (so every computation reads pre-move
    # positions -- the apexes are independent, but staging makes that explicit and order-free).
    spread_cos_min = cos(deg2rad(normal_spread_max_deg))  # skip a pair with dot below this
    moves = Vector{Tuple{Int,SVector{3,Float64}}}()
    for apex in sort!(collect(keys(bases)))               # deterministic order (cosmetic)
        xa = mesh.X[apex]
        triples = bases[apex]

        # Base-face normals oriented away from the apex (= into the material, since the apex
        # is outside), and the mean apex->base edge length h (~ one lattice edge).
        normals = SVector{3,Float64}[]
        h_sum = 0.0
        h_cnt = 0
        for (i, j, k) in triples
            n, ok = face_outward_normal(mesh.X[i], mesh.X[j], mesh.X[k], xa)
            ok || continue
            push!(normals, n)
            h_sum += norm(mesh.X[i] - xa) + norm(mesh.X[j] - xa) + norm(mesh.X[k] - xa)
            h_cnt += 3
        end
        isempty(normals) && continue

        # Fan-spread guard: skip if any two base normals diverge by more than the threshold
        # (their average would not be a meaningful move direction).
        spread_ok = true
        @inbounds for a = 1:(length(normals) - 1), b = (a + 1):length(normals)
            if dot(normals[a], normals[b]) < spread_cos_min
                spread_ok = false
                break
            end
        end
        spread_ok || continue

        n_sum = sum(normals)
        norm(n_sum) < 1e-14 && continue           # opposing normals cancel -> no direction
        n_avg = n_sum / norm(n_sum)
        h = h_sum / h_cnt

        # Move the apex onto phi = 0 along n_avg by bisection of eval_sdf in the [-h, +h]
        # bracket (reused from the relaxation pass). A spike that does not bracket a zero
        # within ~one edge is left as residual.
        p, ok = reproject_to_surface(mesh, xa, n_avg, h, tol)
        ok || continue
        push!(moves, (apex, p))
    end

    # (4) Apply all staged moves at once: snap the apex onto phi = 0 and relabel it on-surface,
    # turning each spike into a quadruple-zero tet the slice will retain.
    for (apex, p) in moves
        mesh.X[apex] = p
        mesh.node_sdf[apex] = 0.0
    end

    @info "Cap recovery: moved $(length(moves)) spike apexes onto phi = 0 " *
          "(of $(length(bases)) +000 apex groups detected)"
    return length(moves)
end
