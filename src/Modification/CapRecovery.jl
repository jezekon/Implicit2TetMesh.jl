# ==============================================================================
# Etapa 11 -- optional convex boundary cap recovery (topology-growing). DEFAULT OFF.
# ==============================================================================
# Where the true surface is convex, the flat boundary triangles of the trimmed mesh
# under-cut it: the mesh boundary chords across the bulge, leaving a thin crescent of
# material un-meshed. This post-pass recovers that material LOCALLY, one boundary tet
# at a time, by a 1:3 split that lifts a new vertex onto the real surface.
#
# A boundary tet ABCD has its face ABC on phi = 0 (the three vertices were warped/cut
# onto the surface) and its apex D strictly inside. We:
#   1. take the face centroid c_f = (A+B+C)/3, and only proceed if it is strictly inside
#      the TRUE field (eval_sdf(c_f) < 0) -- i.e. the flat face really does under-cut a
#      convex bulge (a flat or concave region has eval_sdf(c_f) >= 0 and is left alone);
#   2. project c_f outward (away from D) onto phi = 0 to get P, by the same bracketed
#      bisection the relaxation uses (reproject_to_surface) -- so this is SOURCE-AGNOSTIC
#      (it only reads eval_sdf and works for UnstructuredSDF too);
#   3. split ABCD into the three tets [A,B,P,D], [B,C,P,D], [C,A,P,D] -- the bipyramid on
#      axis D--P. D and P sit on opposite sides of ABC, so the cap adds volume and the
#      boundary now follows the curve.
#
# CONFORMITY IS A CONSTRUCTION, NOT A HOPE. The split lives entirely inside ABCD union the
# cap: the three INTERIOR faces ABD, BCD, CAD are each preserved (still a face of exactly
# one sub-tet, shared with the original neighbour across it), and every edge AB, BC, CA,
# AD, BD, CD survives undivided. The only face that changes is the BOUNDARY face ABC (which
# had no neighbour), replaced by ABP, BCP, CAP. P is a fresh vertex of degree 3 on phi = 0.
# Therefore:
#   - two adjacent boundary tets split INDEPENDENTLY: a shared boundary edge (e.g. AB) is
#     never subdivided by either split, so no hanging node can appear;
#   - watertightness is preserved: each new boundary edge (AP, BP, CP) is shared by exactly
#     two cap faces, and each old boundary edge (AB, BC, CA) keeps its two boundary faces.
#
# SAFETY: every sub-tet is gated by the EXACT orientation predicate is_positively_oriented
# (no tolerance) plus a degenerate-volume floor; if ANY of the three fails, the whole split
# is discarded and the parent tet is kept unchanged. D inside + P outside makes the
# bipyramid genuinely 3D, so well-formed caps pass; the gate only ever rejects, never warps.
#
# DETERMINISTIC: tets are visited in index order, candidates are COLLECTED first and APPLIED
# in a second batch, there is no rand(). Same input -> same output, so caps baselines can be
# frozen in test/.
#
# This is NOT the removed global volume corrector. That moved EXISTING surface nodes onto a
# single SDF level, degrading fidelity and risking spikes/inversions. This INSERTS new nodes
# on the TRUE zero level set, only where the surface is convex, gated on curvature (sagitta)
# and on exact orientation -- a principled, local, conforming refinement of the boundary.

"""
    CapRecoveryOptions(; sagitta_frac = 0.2, bracket_frac = 1.0, max_passes = 1,
                         bisection_tol = 1e-7, min_volume_frac = 0.0)

Configuration for the optional [`recover_boundary_caps!`](@ref) post-pass (Etapa 11).
All defaults are starting points meant to be tuned while measuring.

# Fields
- `sagitta_frac::Float64`: split a boundary face only when the recovered point P lifts off
  the flat face by more than `sagitta_frac * h_face` (`h_face` = mean edge length of the
  face). The sagitta is the GEOMETRIC distance `|P - c_f|`, not `|eval_sdf|` (the latter is
  not a distance for non-distance / unstructured fields). Larger values recover only the
  sharpest bulges; `0.0` recovers every convex under-cut.
- `bracket_frac::Float64`: the outward search half-length is `bracket_frac * h_face`. The
  projection of the centroid onto phi = 0 is sought within `[-h, +h]` along the outward
  normal; a bulge deeper than `h` is left unrecovered (no zero is bracketed -> skipped).
- `max_passes::Int`: how many times to repeat the pass. Each pass exposes the new cap faces
  ABP/BCP/CAP as fresh boundary faces that can themselves be split, refining the bulge
  further. Each cap adds two tets, so passes > 1 grow the element count; default 1.
- `bisection_tol::Float64`: `|eval_sdf|` tolerance for the outward surface projection.
- `min_volume_frac::Float64`: discard a split if any sub-tet's volume is below
  `min_volume_frac * (parent tet volume)`. `0.0` (default) gates only true degeneracy (the
  exact orientation predicate already guarantees a strictly positive volume); raise it to
  refuse caps that would create slivers.
"""
struct CapRecoveryOptions
    sagitta_frac::Float64
    bracket_frac::Float64
    max_passes::Int
    bisection_tol::Float64
    min_volume_frac::Float64

    function CapRecoveryOptions(;
        sagitta_frac::Float64 = 0.2,
        bracket_frac::Float64 = 1.0,
        max_passes::Int = 1,
        bisection_tol::Float64 = 1e-7,
        min_volume_frac::Float64 = 0.0,
    )
        sagitta_frac >= 0.0 || error("sagitta_frac must be >= 0, got $sagitta_frac")
        bracket_frac > 0.0 || error("bracket_frac must be positive, got $bracket_frac")
        max_passes >= 1 || error("max_passes must be >= 1, got $max_passes")
        bisection_tol > 0.0 || error("bisection_tol must be positive, got $bisection_tol")
        min_volume_frac >= 0.0 ||
            error("min_volume_frac must be >= 0, got $min_volume_frac")
        new(sagitta_frac, bracket_frac, max_passes, bisection_tol, min_volume_frac)
    end
end

# ------------------------------------------------------------------------------
# Sub-tet construction with exact orientation + volume gating
# ------------------------------------------------------------------------------

# Sentinel node index standing in for the not-yet-inserted cap apex P inside a sub-tet
# tuple. Real node indices are 1-based, so 0 is unambiguous; Phase C substitutes P's true
# index once it is pushed onto mesh.X.
const CAP_P_SENTINEL = 0

"""
    cap_coord(mesh, P, i) -> SVector{3,Float64}

Coordinate of sub-tet vertex `i`: the cap apex position `P` when `i` is the
[`CAP_P_SENTINEL`](@ref), otherwise the mesh node `mesh.X[i]`.
"""
@inline function cap_coord(mesh::BlockMesh, P::SVector{3,Float64}, i::Int)
    return i == CAP_P_SENTINEL ? P : mesh.X[i]
end

"""
    oriented_subtet(mesh, P, i1, i2, i3, i4, vol_floor) -> Union{Vector{Int},Nothing}

Build one sub-tet `(i1, i2, i3, i4)` (with `P` substituted for the sentinel), GATED:
returns the vertex list reordered to strictly positive orientation, or `nothing` if the
tet is degenerate (volume `<= vol_floor`, or coplanar by the exact predicate). Orientation
is fixed exactly as `remove_inverted_elements!` does -- a single 3<->4 transposition flips
the sign of a non-degenerate tet -- and the sign itself is decided by the shared exact
predicate `is_positively_oriented`.
"""
function oriented_subtet(
    mesh::BlockMesh,
    P::SVector{3,Float64},
    i1::Int,
    i2::Int,
    i3::Int,
    i4::Int,
    vol_floor::Float64,
)
    p1 = cap_coord(mesh, P, i1)
    p2 = cap_coord(mesh, P, i2)
    p3 = cap_coord(mesh, P, i3)
    p4 = cap_coord(mesh, P, i4)
    # Degenerate-volume floor (magnitude only); the SIGN is decided exactly below.
    abs(dot(p2 - p1, cross(p3 - p1, p4 - p1))) / 6.0 > vol_floor || return nothing
    if is_positively_oriented(p1, p2, p3, p4)
        return Int[i1, i2, i3, i4]
    elseif is_positively_oriented(p1, p2, p4, p3)
        return Int[i1, i2, i4, i3]    # one transposition flips a non-degenerate tet
    else
        return nothing               # exactly coplanar: refuse the split
    end
end

"""
    build_cap_subtets(mesh, A, B, C, D, P, vol_floor) -> Union{NTuple{3,Vector{Int}},Nothing}

The three sub-tets of the bipyramid split of `ABCD` about axis `D--P`:
`[A,B,P,D]`, `[B,C,P,D]`, `[C,A,P,D]` (P held as [`CAP_P_SENTINEL`](@ref)), each oriented
positively and volume-gated by [`oriented_subtet`](@ref). Returns `nothing` if ANY sub-tet
fails the gate, so the caller discards the whole split and leaves the parent tet untouched.
"""
function build_cap_subtets(
    mesh::BlockMesh,
    A::Int,
    B::Int,
    C::Int,
    D::Int,
    P::SVector{3,Float64},
    vol_floor::Float64,
)
    s = CAP_P_SENTINEL
    t1 = oriented_subtet(mesh, P, A, B, s, D, vol_floor)
    t1 === nothing && return nothing
    t2 = oriented_subtet(mesh, P, B, C, s, D, vol_floor)
    t2 === nothing && return nothing
    t3 = oriented_subtet(mesh, P, C, A, s, D, vol_floor)
    t3 === nothing && return nothing
    return (t1, t2, t3)
end

# ------------------------------------------------------------------------------
# One pass: collect convex-cap candidates, then apply them as a batch
# ------------------------------------------------------------------------------

"""
    recover_caps_pass!(mesh, opts) -> Int

Run one cap-recovery pass over the current mesh and return the number of caps added (each
turns one tet into three and inserts one surface node). Topology is read from `mesh.IEN`
only -- `mesh.INE` is NOT required and is left stale for the caller to rebuild. See the
file header for the conformity and watertightness guarantees.
"""
function recover_caps_pass!(mesh::BlockMesh, opts::CapRecoveryOptions)
    grid_tol = mesh.grid_tol
    vol_floor = grid_tol * 1e-6                      # same degeneracy floor the pipeline uses

    # --- Phase A: boundary-face incidence map (a face in exactly one tet is on the surface).
    face_count = Dict{NTuple{3,Int},Int}()
    @inbounds for tet in mesh.IEN
        a, b, c, d = tet[1], tet[2], tet[3], tet[4]
        for (x, y, z) in ((a, b, c), (a, b, d), (a, c, d), (b, c, d))
            k = face_key(x, y, z)
            face_count[k] = get(face_count, k, 0) + 1
        end
    end

    # --- Phase B: collect candidates (no mutation). Deterministic index order.
    cand_tet = Int[]                                 # parent tet indices (ascending)
    cand_P = SVector{3,Float64}[]                    # recovered apex positions
    cand_sub = NTuple{3,Vector{Int}}[]               # the three oriented sub-tets each
    @inbounds for (ti, tet) in enumerate(mesh.IEN)
        a, b, c, d = tet[1], tet[2], tet[3], tet[4]

        # (1) exactly one boundary face ABC; its opposite vertex is the interior apex D.
        nbound = 0
        A = B = C = D = 0
        for (x, y, z, o) in ((a, b, c, d), (a, b, d, c), (a, c, d, b), (b, c, d, a))
            if face_count[face_key(x, y, z)] == 1
                nbound += 1
                A, B, C, D = x, y, z, o
            end
        end
        nbound == 1 || continue

        # (2) the boundary face sits on phi = 0; the apex is strictly inside.
        (abs(mesh.node_sdf[A]) <= grid_tol &&
         abs(mesh.node_sdf[B]) <= grid_tol &&
         abs(mesh.node_sdf[C]) <= grid_tol) || continue
        mesh.node_sdf[D] < -grid_tol || continue

        pA = mesh.X[A]; pB = mesh.X[B]; pC = mesh.X[C]; pD = mesh.X[D]

        # (3) convex under-cut only: the face centroid must lie strictly INSIDE the true
        #     field. A flat/concave region has eval_sdf(c_f) >= 0 and is left alone.
        c_f = (pA + pB + pC) / 3.0
        eval_sdf(mesh, c_f) < -grid_tol || continue

        # (4) outward face normal (away from the interior apex D).
        n_out, okn = face_outward_normal(pA, pB, pC, pD)
        okn || continue

        # (5) project the centroid onto phi = 0 along the outward normal.
        h_face = (norm(pB - pA) + norm(pC - pB) + norm(pA - pC)) / 3.0
        P, okp = reproject_to_surface(
            mesh, c_f, n_out, opts.bracket_frac * h_face, opts.bisection_tol)
        okp || continue

        # (6) skip negligible bulges -- GEOMETRIC sagitta (valid for non-distance fields).
        norm(P - c_f) >= opts.sagitta_frac * h_face || continue

        # (7) build + gate the three sub-tets; reject the whole split if any is degenerate.
        parent_vol = abs(dot(pB - pA, cross(pC - pA, pD - pA))) / 6.0
        subs = build_cap_subtets(mesh, A, B, C, D, P, max(vol_floor, opts.min_volume_frac * parent_vol))
        subs === nothing && continue

        push!(cand_tet, ti)
        push!(cand_P, P)
        push!(cand_sub, subs)
    end

    isempty(cand_tet) && return 0

    # --- Phase C: apply the batch deterministically.
    ncand = length(cand_tet)
    catpos = zeros(Int, length(mesh.IEN))            # tet index -> candidate slot (0 = none)
    for i = 1:ncand
        catpos[cand_tet[i]] = i
    end

    # Insert the recovered apices; each lies exactly on phi = 0 (node_sdf = 0).
    pidx = Vector{Int}(undef, ncand)
    for i = 1:ncand
        push!(mesh.X, cand_P[i])
        push!(mesh.node_sdf, 0.0)
        pidx[i] = length(mesh.X)
    end

    # Replace each marked tet by its three sub-tets; copy the rest. Order preserved.
    new_IEN = Vector{Vector{Int}}()
    sizehint!(new_IEN, length(mesh.IEN) + 2 * ncand)
    @inbounds for (ti, tet) in enumerate(mesh.IEN)
        p = catpos[ti]
        if p == 0
            push!(new_IEN, tet)
        else
            pi = pidx[p]
            for sub in cand_sub[p]
                push!(new_IEN, Int[
                    sub[1] == CAP_P_SENTINEL ? pi : sub[1],
                    sub[2] == CAP_P_SENTINEL ? pi : sub[2],
                    sub[3] == CAP_P_SENTINEL ? pi : sub[3],
                    sub[4] == CAP_P_SENTINEL ? pi : sub[4],
                ])
            end
        end
    end
    mesh.IEN = new_IEN
    return ncand
end

# ------------------------------------------------------------------------------
# Driver
# ------------------------------------------------------------------------------

"""
    recover_boundary_caps!(mesh::BlockMesh, opts::CapRecoveryOptions) -> mesh

Optional convex boundary cap recovery (Etapa 11, default OFF). Run on the final mesh, AFTER
the last `update_connectivity!` (so `mesh.node_sdf` and the boundary are valid) and BEFORE
relaxation / export. Where the flat boundary under-cuts a convex bulge, a boundary tet is
split 1:3 by inserting a vertex projected onto the true zero level set, recovering the
missing material while preserving conformity, watertightness, and positive orientation by
construction (see the file header).

Like `remove_inverted_elements!`, this changes topology and does NOT rebuild the inverse
connectivity `mesh.INE` or compact nodes -- the caller refreshes connectivity right after
(`update_connectivity!`). It only reads the field through `eval_sdf`, so it works
identically for structured and unstructured sources. Returns the modified mesh.
"""
function recover_boundary_caps!(mesh::BlockMesh, opts::CapRecoveryOptions)
    isempty(mesh.IEN) && return mesh
    @info "Recovering convex boundary caps " *
          "(sagitta_frac = $(opts.sagitta_frac), bracket_frac = $(opts.bracket_frac), " *
          "max_passes = $(opts.max_passes))..."

    total = 0
    for pass = 1:opts.max_passes
        added = recover_caps_pass!(mesh, opts)
        total += added
        @info "  pass $pass: recovered $added caps (+$(2 * added) tets, +$added nodes)"
        added == 0 && break
    end
    @info "  Total: $total caps recovered (+$(2 * total) tets, +$total nodes)."
    return mesh
end
