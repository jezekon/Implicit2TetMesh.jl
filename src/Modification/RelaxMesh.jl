# ==============================================================================
# Etapa 10 -- optional, gated vertex relaxation (fixed topology). DEFAULT OFF.
# ==============================================================================
# The A15 interior is already uniform and near-optimal; all element-size variance
# and all bad dihedral angles live in the 1-2 element layers created by warp + trim.
# This post-pass improves THAT band without ever changing the topology (no flips,
# collapses, or insertions): it only MOVES vertices. Two objectives share one
# skeleton -- PROPOSE a new vertex position, GATE it on local element quality,
# PROJECT surface vertices back to phi = 0:
#
#   mode = :uniform  -- DistMesh size-equalizing springs (Persson & Strang 2004):
#                       FE users want similar-sized elements near the boundary.
#   mode = :quality  -- quartet-style maximin smoothing (optimize_tet_mesh.cpp):
#                       actively lift the worst dihedral angles.
#
# SAFETY IS A CONSTRUCTION, NOT A HOPE. The current position always competes, and a
# move is accepted only if it does not lower the local minimum quality. Therefore the
# global minimum quality is monotonically non-decreasing, so "min dihedral >= pre-pass"
# is a hard, testable guarantee. Surface vertices stay EXACTLY on phi = 0 at all times
# (tangential proposal + re-projection along the angle-weighted pseudo-normal), so this
# is NOT the removed volume corrector -- it never moves a node off the zero level; the
# volume changes only by chord redistribution.
#
# DETERMINISTIC: no rand() (quartet uses rand; we do not). Fixed-order Gauss-Seidel
# sweeps over a sorted active set, fixed candidate patterns. Same input -> same output,
# so relaxed baselines can be frozen in test/.

"""
    RelaxOptions(; mode = :uniform, max_sweeps = 10, band = 2, omega = 0.3,
                   sizing = :uniform, alpha = 0.5, h_min = 0.0, h_max = 0.0,
                   grading = 0.3, frozen = Int[], bisection_tol = 1e-7)

Configuration for the optional [`relax_mesh!`](@ref) post-pass. All defaults are
starting points meant to be tuned while measuring.

# Fields
- `mode::Symbol`: `:uniform` (size-equalizing springs) or `:quality` (maximin
  dihedral smoothing).
- `max_sweeps::Int`: maximum number of Gauss-Seidel sweeps over the active set.
- `band::Int`: the active set is the boundary-surface vertices plus `band` topological
  rings inward (the A15 interior is already optimal). Use a very large value
  (e.g. `typemax(Int)`) to relax every vertex -- an experiment-only override.
- `omega::Float64`: spring step factor for `mode = :uniform`.
- `sizing::Symbol`: `:uniform` (one global target edge length) or `:curvature`
  (per-vertex target concentrating elements where the surface bends). Only meaningful
  for `mode = :uniform`.
- `alpha::Float64`: curvature sizing scale; target length ~ `alpha / curvature`.
- `h_min::Float64`, `h_max::Float64`: clamps on the curvature target length. `0.0`
  means "auto" (`0.5 * grid_step` and `4.0 * grid_step` respectively).
- `grading::Float64`: gradient-limit bound `g` for curvature sizing (Persson 2006):
  neighbouring target lengths may differ by at most `g * edge_length`.
- `frozen::Vector{Int}`: externally supplied freeze-set (node indices never moved);
  the Etapa-9 functional-surface hook plugs in here.
- `bisection_tol::Float64`: `|eval_sdf|` tolerance for the surface re-projection.
"""
struct RelaxOptions
    mode::Symbol
    max_sweeps::Int
    band::Int
    omega::Float64
    sizing::Symbol
    alpha::Float64
    h_min::Float64
    h_max::Float64
    grading::Float64
    frozen::Vector{Int}
    bisection_tol::Float64

    function RelaxOptions(;
        mode::Symbol = :uniform,
        max_sweeps::Int = 10,
        band::Int = 2,
        omega::Float64 = 0.3,
        sizing::Symbol = :uniform,
        alpha::Float64 = 0.5,
        h_min::Float64 = 0.0,
        h_max::Float64 = 0.0,
        grading::Float64 = 0.3,
        frozen::AbstractVector{<:Integer} = Int[],
        bisection_tol::Float64 = 1e-7,
    )
        mode === :uniform || mode === :quality ||
            error("Invalid relax mode: $mode. Use :uniform or :quality.")
        sizing === :uniform || sizing === :curvature ||
            error("Invalid relax sizing: $sizing. Use :uniform or :curvature.")
        max_sweeps >= 1 || error("max_sweeps must be >= 1, got $max_sweeps")
        band >= 0 || error("band must be >= 0, got $band")
        omega > 0.0 || error("omega must be positive, got $omega")
        alpha > 0.0 || error("alpha must be positive, got $alpha")
        h_min >= 0.0 || error("h_min must be >= 0, got $h_min")
        h_max >= 0.0 || error("h_max must be >= 0, got $h_max")
        grading > 0.0 || error("grading must be positive, got $grading")
        bisection_tol > 0.0 || error("bisection_tol must be positive, got $bisection_tol")
        new(mode, max_sweeps, band, omega, sizing, alpha, h_min, h_max, grading,
            Vector{Int}(frozen), bisection_tol)
    end
end

# ------------------------------------------------------------------------------
# Dihedral-angle quality metric (cosine extremes of the six dihedral angles)
# ------------------------------------------------------------------------------

"""
    face_outward_normal(pa, pb, pc, popp) -> (n, ok)

Unit normal of triangle `(pa, pb, pc)` oriented OUTWARD, i.e. away from the opposite
tetrahedron vertex `popp`. `ok` is false for a degenerate (zero-area) face.
"""
function face_outward_normal(
    pa::SVector{3,Float64},
    pb::SVector{3,Float64},
    pc::SVector{3,Float64},
    popp::SVector{3,Float64},
)
    n = cross(pb - pa, pc - pa)
    l = norm(n)
    l < 1e-14 && return (zero(SVector{3,Float64}), false)
    n /= l
    if dot(n, popp - (pa + pb + pc) / 3.0) > 0.0
        n = -n
    end
    return (n, true)
end

"""
    tet_dihedral_cos_extremes(p1, p2, p3, p4) -> (cmin, cmax, sin2min)

Cosine extremes and minimum squared-sine of the six interior dihedral angles of the
tetrahedron, the quality scalars the relaxation gate works with. With OUTWARD face
normals the cosine of the dihedral on a shared edge is `-dot(n_a, n_b)`, so a sharper
angle (toward 0 deg) has cosine toward +1 and a wider angle (toward 180 deg) has cosine
toward -1. Therefore:
  - `cmax` = the largest cosine = the SHARPEST dihedral (the sliver tail);
  - `cmin` = the smallest cosine = the WIDEST dihedral (the cap tail);
  - `sin2min` = `min(1 - cos^2)` = quartet's min-sine-of-dihedral quality (squared).
A two-sided gate that keeps `cmax` from rising and `cmin` from falling protects BOTH
tails -- guaranteeing min dihedral does not drop AND max dihedral does not grow. The
six edges correspond one-to-one to the six face pairs. A degenerate tet returns the
worst values `(-1, 1, 0)`; orientation is decided separately and exactly by
`is_positively_oriented`.
"""
function tet_dihedral_cos_extremes(
    p1::SVector{3,Float64},
    p2::SVector{3,Float64},
    p3::SVector{3,Float64},
    p4::SVector{3,Float64},
)
    nA, okA = face_outward_normal(p2, p3, p4, p1)   # opposite p1
    nB, okB = face_outward_normal(p1, p3, p4, p2)   # opposite p2
    nC, okC = face_outward_normal(p1, p2, p4, p3)   # opposite p3
    nD, okD = face_outward_normal(p1, p2, p3, p4)   # opposite p4
    (okA && okB && okC && okD) || return (-1.0, 1.0, 0.0)

    # Cosine of the interior dihedral on each of the six edges (the two faces sharing it).
    cmin = Inf
    cmax = -Inf
    sin2min = Inf
    for (na, nb) in ((nC, nD), (nB, nD), (nB, nC), (nA, nD), (nA, nC), (nA, nB))
        c = clamp(-dot(na, nb), -1.0, 1.0)
        cmin = min(cmin, c)
        cmax = max(cmax, c)
        sin2min = min(sin2min, 1.0 - c * c)
    end
    return (cmin, cmax, sin2min)
end

# ------------------------------------------------------------------------------
# Surface re-projection onto phi = 0 (bisection along the pseudo-normal)
# ------------------------------------------------------------------------------

"""
    reproject_to_surface(mesh, x, n, h, tol) -> (point, converged)

Move `x` onto the zero level set `phi = 0` along the unit pseudo-normal `n`, by
bisection of `eval_sdf` inside the clamped bracket `[-h, +h]`. Clamping the bracket to
half the local edge length keeps the projection from jumping to the opposite sheet of a
thin wall. The search direction is the angle-weighted pseudo-normal -- NOT the SDF
gradient (which was removed in audit Part B and is invalid for non-distance fields).

Returns `(projected_point, true)` when a zero is bracketed and refined to within `tol`,
and `(x, false)` when no sign change is found in the bracket; the caller then treats the
candidate as invalid so a surface vertex is never left off the surface.
"""
function reproject_to_surface(
    mesh::BlockMesh,
    x::SVector{3,Float64},
    n::SVector{3,Float64},
    h::Float64,
    tol::Float64,
)
    f0 = eval_sdf(mesh, x)
    abs(f0) <= tol && return (x, true)

    fa = eval_sdf(mesh, x - h * n)   # value at the -h end
    fb = eval_sdf(mesh, x + h * n)   # value at the +h end

    # Pick a sub-bracket that contains a sign change, preferring the side nearest x.
    local lo::Float64, hi::Float64, flo::Float64
    if (f0 < 0) != (fb < 0)
        lo = 0.0;  hi = h;   flo = f0
    elseif (f0 < 0) != (fa < 0)
        lo = -h;   hi = 0.0; flo = fa
    elseif (fa < 0) != (fb < 0)
        lo = -h;   hi = h;   flo = fa
    else
        return (x, false)            # no zero within the clamped bracket
    end

    # Bisection. 60 halvings of a bracket ~ edge length drive |t| well below any tol.
    for _ = 1:60
        mid = 0.5 * (lo + hi)
        fmid = eval_sdf(mesh, x + mid * n)
        if abs(fmid) <= tol
            return (x + mid * n, true)
        end
        if (flo < 0) != (fmid < 0)
            hi = mid
        else
            lo = mid
            flo = fmid
        end
    end
    mid = 0.5 * (lo + hi)
    return (x + mid * n, true)
end

# ------------------------------------------------------------------------------
# Shared infrastructure, built once per relax_mesh! call
# ------------------------------------------------------------------------------

"""
    RelaxInfra

Per-call infrastructure for the relaxation sweeps. Built once because the topology
never changes during the pass, so `mesh.INE`, `mesh.IEN` and the boundary
triangulation stay valid throughout -- only vertex positions move.

# Fields (all indexed by global node id, 1..length(mesh.X))
- `is_surface::BitVector`: vertex lies on the boundary surface.
- `is_movable::BitVector`: vertex is in the active band, not on a non-manifold pinch,
  and not in the user freeze-set.
- `movable::Vector{Int}`: the movable node ids in ascending order (sweep order).
- `pnormal::Vector{SVector{3,Float64}}`: angle-weighted outward pseudo-normal at each
  surface vertex (zero elsewhere).
- `local_h::Vector{Float64}`: local edge-length scale at each surface vertex (mean of
  its boundary-edge lengths); half of it is the re-projection bracket.
- `L0::Vector{Float64}`: per-vertex target edge length for the springs.
- `interior_nb::Dict{Int,Vector{Int}}`: tet-edge neighbours of each movable INTERIOR
  vertex (the full spring stencil).
- `surface_nb::Dict{Int,Vector{Int}}`: boundary-edge neighbours of each movable SURFACE
  vertex (the surface-only spring stencil).
"""
struct RelaxInfra
    is_surface::BitVector
    is_movable::BitVector
    movable::Vector{Int}
    pnormal::Vector{SVector{3,Float64}}
    local_h::Vector{Float64}
    L0::Vector{Float64}
    interior_nb::Dict{Int,Vector{Int}}
    surface_nb::Dict{Int,Vector{Int}}
end

"""
    vertex_neighbors(mesh, v) -> Vector{Int}

Distinct tet-edge neighbours of vertex `v`, read from the inverse connectivity
`mesh.INE` (the union of the other three vertices of every incident tet). Used both for
the band BFS and for the interior spring stencil.
"""
function vertex_neighbors(mesh::BlockMesh, v::Int)
    nb = Set{Int}()
    @inbounds for e in mesh.INE[v]
        for w in mesh.IEN[e]
            w != v && push!(nb, w)
        end
    end
    return sort!(collect(nb))
end

"""
    build_relax_infra(mesh, opts) -> RelaxInfra

Assemble the boundary triangulation, pseudo-normals, freeze-set, active band and spring
stencils for one relaxation pass. See [`RelaxInfra`](@ref) for the contents.
"""
function build_relax_infra(mesh::BlockMesh, opts::RelaxOptions)
    n_nodes = length(mesh.X)

    # --- Boundary faces (incidence 1), remembering one outward-oriented triangle each.
    # A face shared by two tets is interior; a face in exactly one tet is on the surface.
    # For each face we keep an example (a, b, c, opposite) so we can orient it outward.
    face_count = Dict{NTuple{3,Int},Int}()
    face_repr = Dict{NTuple{3,Int},NTuple{4,Int}}()
    @inbounds for tet in mesh.IEN
        a, b, c, d = tet[1], tet[2], tet[3], tet[4]
        for (x, y, z, o) in ((a, b, c, d), (a, b, d, c), (a, c, d, b), (b, c, d, a))
            key = face_key(x, y, z)
            face_count[key] = get(face_count, key, 0) + 1
            face_repr[key] = (x, y, z, o)
        end
    end

    is_surface = falses(n_nodes)
    pnormal = fill(zero(SVector{3,Float64}), n_nodes)
    edge_len_sum = zeros(Float64, n_nodes)
    edge_len_cnt = zeros(Int, n_nodes)
    surface_nb_set = Dict{Int,Set{Int}}()
    boundary_edge_count = Dict{NTuple{2,Int},Int}()

    for (key, cnt) in face_count
        cnt == 1 || continue
        (x, y, z, o) = face_repr[key]
        # Orient the triangle so its normal points away from the opposite vertex o.
        p1 = mesh.X[x]; p2 = mesh.X[y]; p3 = mesh.X[z]
        nf = cross(p2 - p1, p3 - p1)
        nlen = norm(nf)
        nlen < 1e-14 && continue
        nf /= nlen
        if dot(nf, mesh.X[o] - (p1 + p2 + p3) / 3.0) > 0.0
            nf = -nf   # flip to outward
        end

        tri = (x, y, z)
        # Angle-weighted pseudo-normal contribution + boundary edges + edge lengths.
        for (vi, va, vb) in ((x, y, z), (y, x, z), (z, x, y))
            is_surface[vi] = true
            e1 = mesh.X[va] - mesh.X[vi]
            e2 = mesh.X[vb] - mesh.X[vi]
            l1 = norm(e1); l2 = norm(e2)
            if l1 > 1e-14 && l2 > 1e-14
                ang = acos(clamp(dot(e1, e2) / (l1 * l2), -1.0, 1.0))
                pnormal[vi] += ang * nf
            end
            nbset = get!(surface_nb_set, vi, Set{Int}())
            push!(nbset, va)
            push!(nbset, vb)
            edge_len_sum[vi] += l1
            edge_len_cnt[vi] += 1
        end
        for (u, w) in ((x, y), (x, z), (y, z))
            ek = (min(u, w), max(u, w))
            boundary_edge_count[ek] = get(boundary_edge_count, ek, 0) + 1
        end
    end

    # Normalise the pseudo-normals; record the local edge scale.
    local_h = zeros(Float64, n_nodes)
    for v = 1:n_nodes
        is_surface[v] || continue
        pl = norm(pnormal[v])
        pnormal[v] = pl > 1e-14 ? pnormal[v] / pl : zero(SVector{3,Float64})
        local_h[v] = edge_len_cnt[v] > 0 ? edge_len_sum[v] / edge_len_cnt[v] : mesh.grid_step
    end

    # --- Freeze pinch vertices: a manifold boundary edge is shared by exactly two
    # boundary triangles. Vertices on edges with any other count sit on a non-manifold
    # pinch or a crack, where the pseudo-normal is meaningless -- never move them.
    pinch = falses(n_nodes)
    for (ek, c) in boundary_edge_count
        if c != 2
            pinch[ek[1]] = true
            pinch[ek[2]] = true
        end
    end

    # --- Active band: surface vertices (ring 0) plus `band` rings inward via BFS over
    # tet-edge neighbours. A huge `band` relaxes the whole mesh.
    in_band = falses(n_nodes)
    frontier = Int[]
    for v = 1:n_nodes
        if is_surface[v]
            in_band[v] = true
            push!(frontier, v)
        end
    end
    for _ = 1:opts.band
        isempty(frontier) && break
        next_frontier = Int[]
        for v in frontier
            for w in vertex_neighbors(mesh, v)
                if !in_band[w]
                    in_band[w] = true
                    push!(next_frontier, w)
                end
            end
        end
        frontier = next_frontier
    end

    # --- Freeze-set from options, then assemble the movable list and spring stencils.
    user_frozen = falses(n_nodes)
    for v in opts.frozen
        1 <= v <= n_nodes && (user_frozen[v] = true)
    end

    is_movable = falses(n_nodes)
    interior_nb = Dict{Int,Vector{Int}}()
    surface_nb = Dict{Int,Vector{Int}}()
    movable = Int[]
    for v = 1:n_nodes
        (in_band[v] && !pinch[v] && !user_frozen[v]) || continue
        is_movable[v] = true
        push!(movable, v)
        if is_surface[v]
            surface_nb[v] = sort!(collect(surface_nb_set[v]))
        else
            interior_nb[v] = vertex_neighbors(mesh, v)
        end
    end

    # --- Target edge length L0 (springs). :uniform -> one global median over active-band
    # edges; :curvature -> per-vertex (filled in below).
    L0 = fill(default_target_length(mesh, in_band), n_nodes)
    if opts.mode === :uniform && opts.sizing === :curvature
        fill_curvature_sizing!(L0, mesh, opts, is_surface, is_movable, pnormal,
                               local_h, surface_nb, interior_nb)
    end

    return RelaxInfra(is_surface, is_movable, movable, pnormal, local_h, L0,
                      interior_nb, surface_nb)
end

"""
    default_target_length(mesh, in_band) -> Float64

Median edge length over the active band -- the global target edge length `L0` for the
uniform springs. Falls back to `grid_step` if the band has no edges.
"""
function default_target_length(mesh::BlockMesh, in_band::BitVector)
    lengths = Float64[]
    seen = Set{NTuple{2,Int}}()
    @inbounds for tet in mesh.IEN
        for (a, b) in ((tet[1], tet[2]), (tet[1], tet[3]), (tet[1], tet[4]),
                       (tet[2], tet[3]), (tet[2], tet[4]), (tet[3], tet[4]))
            (in_band[a] || in_band[b]) || continue
            ek = (min(a, b), max(a, b))
            if !(ek in seen)
                push!(seen, ek)
                push!(lengths, norm(mesh.X[a] - mesh.X[b]))
            end
        end
    end
    isempty(lengths) && return mesh.grid_step
    return median(lengths)
end

# ------------------------------------------------------------------------------
# Curvature sizing (optional, :uniform mode only) -- Persson 2006 gradient limiting
# ------------------------------------------------------------------------------

"""
    fill_curvature_sizing!(L0, mesh, opts, is_surface, is_movable, pnormal, local_h,
                           surface_nb, interior_nb)

Fill a per-vertex target edge length that CONCENTRATES elements where the surface bends.
For each surface vertex the discrete curvature `kappa` is the largest angle between its
pseudo-normal and a surface-neighbour's pseudo-normal, divided by that edge length (no
SDF derivatives); the raw target is `clamp(alpha / kappa, h_min, h_max)`. The field is
then gradient-limited over the active-band edge graph so neighbouring targets differ by
at most `g * edge_length` (discrete `|grad h| <= g`, Persson 2006), and the limited
length is propagated to interior band vertices.

Honest expectation: with fixed topology this yields only MILD concentration (edge ratios
~1.5-2x), not true refinement -- real refinement is Etapa 6's octree.
"""
function fill_curvature_sizing!(
    L0::Vector{Float64},
    mesh::BlockMesh,
    opts::RelaxOptions,
    is_surface::BitVector,
    is_movable::BitVector,
    pnormal::Vector{SVector{3,Float64}},
    local_h::Vector{Float64},
    surface_nb::Dict{Int,Vector{Int}},
    interior_nb::Dict{Int,Vector{Int}},
)
    h_min = opts.h_min > 0.0 ? opts.h_min : 0.5 * mesh.grid_step
    h_max = opts.h_max > 0.0 ? opts.h_max : 4.0 * mesh.grid_step

    # Raw curvature target at movable surface vertices.
    for v in keys(surface_nb)
        nv = pnormal[v]
        norm(nv) < 1e-14 && (L0[v] = clamp(local_h[v], h_min, h_max); continue)
        kappa = 0.0
        for j in surface_nb[v]
            nj = pnormal[j]
            norm(nj) < 1e-14 && continue
            elen = norm(mesh.X[j] - mesh.X[v])
            elen < 1e-14 && continue
            ang = acos(clamp(dot(nv, nj), -1.0, 1.0))
            kappa = max(kappa, ang / elen)
        end
        L0[v] = kappa > 1e-14 ? clamp(opts.alpha / kappa, h_min, h_max) : h_max
    end

    # Gradient limiting: iterate h_i <- min(h_i, h_j + g*|e_ij|) over band edges until
    # the field stops changing (a few passes -- it is a shortest-path-like relaxation).
    g = opts.grading
    for _ = 1:50
        changed = false
        for v in keys(surface_nb)
            for j in surface_nb[v]
                is_movable[j] || continue
                elen = norm(mesh.X[j] - mesh.X[v])
                bound = L0[j] + g * elen
                if L0[v] > bound
                    L0[v] = bound
                    changed = true
                end
            end
        end
        changed || break
    end

    # Propagate the limited length inward to interior band vertices (mean of movable
    # neighbours' targets), so the springs there match the surface sizing.
    for v in keys(interior_nb)
        s = 0.0; c = 0
        for j in interior_nb[v]
            if is_movable[j]
                s += L0[j]; c += 1
            end
        end
        c > 0 && (L0[v] = s / c)
    end
    return L0
end

# ------------------------------------------------------------------------------
# The gate
# ------------------------------------------------------------------------------

"""
    local_quality(mesh, v, elems) -> (valid, cmin_star, cmax_star, sin2_star)

The quality of vertex `v`'s star, read at the CURRENT `mesh.X`. `valid` is false the
moment any incident tet fails the EXACT orientation predicate `is_positively_oriented`
(no tolerance). Otherwise, over all incident tets:
  - `cmax_star` is the largest dihedral cosine (the SHARPEST angle in the star);
  - `cmin_star` is the smallest dihedral cosine (the WIDEST angle in the star);
  - `sin2_star` is the smallest squared sine (quartet's min-sine quality of the star).
"""
function local_quality(mesh::BlockMesh, v::Int, elems::Vector{Int})
    cmin_star = Inf
    cmax_star = -Inf
    sin2_star = Inf
    @inbounds for e in elems
        tet = mesh.IEN[e]
        is_positively_oriented(mesh, tet) || return (false, -1.0, 1.0, 0.0)
        (cmin, cmax, sin2) = tet_dihedral_cos_extremes(
            mesh.X[tet[1]], mesh.X[tet[2]], mesh.X[tet[3]], mesh.X[tet[4]])
        cmin_star = min(cmin_star, cmin)
        cmax_star = max(cmax_star, cmax)
        sin2_star = min(sin2_star, sin2)
    end
    return (true, cmin_star, cmax_star, sin2_star)
end

"""
    eval_candidate(mesh, v, xnew, elems, is_interior) -> (valid, cmin, cmax, sin2)

Score the candidate position `xnew` for vertex `v` WITHOUT committing: temporarily place
`v` at `xnew`, measure its star quality, then restore. `valid` requires (a) every
incident tet positively oriented (exact) and (b) for an INTERIOR vertex a strictly
inside position (`eval_sdf < 0`), keeping the node-SDF invariant meaningful and
protecting thin features. The caller additionally applies the two-sided dihedral gate.
"""
function eval_candidate(
    mesh::BlockMesh,
    v::Int,
    xnew::SVector{3,Float64},
    elems::Vector{Int},
    is_interior::Bool,
)
    old = mesh.X[v]
    mesh.X[v] = xnew
    valid, cmin, cmax, sin2 = local_quality(mesh, v, elems)
    if valid && is_interior && eval_sdf(mesh, xnew) >= 0.0
        valid = false
    end
    mesh.X[v] = old
    return (valid, cmin, cmax, sin2)
end

"""
    quality_not_worse(cmin_cand, cmax_cand, cmin_cur, cmax_cur) -> Bool

The two-sided dihedral gate: the candidate is accepted only if its star's sharpest angle
is no sharper (`cmax_cand <= cmax_cur`) AND its widest angle is no wider
(`cmin_cand >= cmin_cur`). This is what makes "min dihedral does not drop" and "max
dihedral does not grow" hard guarantees: an accepted move never worsens either tail, so
the global extremes are monotone across the whole pass.
"""
function quality_not_worse(cmin_cand, cmax_cand, cmin_cur, cmax_cur)
    return cmax_cand <= cmax_cur && cmin_cand >= cmin_cur
end

# ------------------------------------------------------------------------------
# Per-vertex relaxation (one Gauss-Seidel update)
# ------------------------------------------------------------------------------

"""
    spring_proposal(mesh, v, neighbors, L0v, omega) -> SVector{3,Float64}

DistMesh size-equalizing displacement at `v`: `omega * sum_j (x_j - x_v)(1 - L0v/|e|)`
over the given neighbour stencil. Short edges push `v` away (compression), long edges
pull it in (tension); the equilibrium is uniform edge lengths around `v`.
"""
function spring_proposal(
    mesh::BlockMesh,
    v::Int,
    neighbors::Vector{Int},
    L0v::Float64,
    omega::Float64,
)
    xv = mesh.X[v]
    d = zero(SVector{3,Float64})
    @inbounds for j in neighbors
        e = mesh.X[j] - xv
        L = norm(e)
        L < 1e-14 && continue
        d += e * (1.0 - L0v / L)
    end
    return omega * d
end

"""
    relax_vertex_uniform!(mesh, v, infra, opts, elems) -> (moved, disp)

One uniform-mode (spring) update of vertex `v`. Interior vertices take the raw spring
displacement; surface vertices project the displacement into the tangent plane and
re-project onto `phi = 0` along the pseudo-normal. The move is gated (orientation +
non-decreasing local quality + interior-inside) and back-tracked over step factors
`{1, 1/2, 1/4}`; if none is accepted the vertex stays put. Returns whether it moved and
the distance moved.
"""
function relax_vertex_uniform!(
    mesh::BlockMesh,
    v::Int,
    infra::RelaxInfra,
    opts::RelaxOptions,
    elems::Vector{Int},
)
    old = mesh.X[v]
    _, cmin_cur, cmax_cur, _ = local_quality(mesh, v, elems)
    is_surf = infra.is_surface[v]

    if is_surf
        n = infra.pnormal[v]
        norm(n) < 1e-14 && return (false, 0.0)   # no usable normal -> skip
        d = spring_proposal(mesh, v, infra.surface_nb[v], infra.L0[v], opts.omega)
        d_tan = d - dot(d, n) * n                  # tangential component only
        bracket = 0.5 * infra.local_h[v]
        for s in (1.0, 0.5, 0.25)
            x_try = old + s * d_tan
            x_proj, ok = reproject_to_surface(mesh, x_try, n, bracket, opts.bisection_tol)
            ok || continue
            valid, cmin, cmax, _ = eval_candidate(mesh, v, x_proj, elems, false)
            if valid && quality_not_worse(cmin, cmax, cmin_cur, cmax_cur)
                mesh.X[v] = x_proj          # surface node stays on phi = 0
                return (true, norm(x_proj - old))
            end
        end
        return (false, 0.0)
    else
        d = spring_proposal(mesh, v, infra.interior_nb[v], infra.L0[v], opts.omega)
        for s in (1.0, 0.5, 0.25)
            x_try = old + s * d
            valid, cmin, cmax, _ = eval_candidate(mesh, v, x_try, elems, true)
            if valid && quality_not_worse(cmin, cmax, cmin_cur, cmax_cur)
                mesh.X[v] = x_try
                mesh.node_sdf[v] = eval_sdf(mesh, x_try)   # refresh moved interior node
                return (true, norm(x_try - old))
            end
        end
        return (false, 0.0)
    end
end

"""
    relax_vertex_quality!(mesh, v, infra, opts, elems) -> (moved, qmin)

One quality-mode (maximin) update of vertex `v`. Scores a FIXED candidate pattern around
the current position with quartet's metric (here `tet_min_sin2_dihedral`) and moves to
the argmax, breaking ties by the smallest displacement. The current position always
competes, so the local minimum quality never decreases. Interior patterns are cube
corners + axis points at radius `P` and `P/2`; surface patterns are the analogous 2D
pattern in the tangent plane, each re-projected to `phi = 0`. Returns whether it moved
and the local minimum quality at the chosen position.
"""
function relax_vertex_quality!(
    mesh::BlockMesh,
    v::Int,
    infra::RelaxInfra,
    opts::RelaxOptions,
    elems::Vector{Int},
)
    old = mesh.X[v]
    _, cmin_cur, cmax_cur, sin2_cur = local_quality(mesh, v, elems)
    P = 0.5 * mesh.grid_step
    is_surf = infra.is_surface[v]

    # Among candidates that pass the two-sided gate (neither tail worse than now),
    # maximise the min-sine quality; tie -> smallest displacement. The current position
    # always competes, so a move is taken only when it strictly improves min-sine.
    best_sin2 = sin2_cur
    best_x = old
    best_d = 0.0

    if is_surf
        n = infra.pnormal[v]
        norm(n) < 1e-14 && return (false, sin2_cur)
        # Tangent basis (quartet's choice), then 2D corner + axis pattern at P and P/2.
        u = abs(n[1]) > 0.5 ? SVector(-n[2], n[1], 0.0) : SVector(0.0, -n[3], n[2])
        u = u / norm(u)
        w = cross(n, u)
        bracket = 0.5 * infra.local_h[v]
        for r in (P, 0.5 * P)
            for (cu, cw) in ((r, r), (-r, r), (r, -r), (-r, -r),
                             (r, 0.0), (-r, 0.0), (0.0, r), (0.0, -r))
                x_try = old + cu * u + cw * w
                x_proj, ok = reproject_to_surface(mesh, x_try, n, bracket, opts.bisection_tol)
                ok || continue
                valid, cmin, cmax, sin2 = eval_candidate(mesh, v, x_proj, elems, false)
                (valid && quality_not_worse(cmin, cmax, cmin_cur, cmax_cur)) || continue
                dmove = norm(x_proj - old)
                if sin2 > best_sin2 || (sin2 == best_sin2 && dmove < best_d)
                    best_sin2 = sin2; best_x = x_proj; best_d = dmove
                end
            end
        end
        if best_d > 0.0
            mesh.X[v] = best_x                 # stays on phi = 0
            return (true, best_sin2)
        end
        return (false, sin2_cur)
    else
        for r in (P, 0.5 * P)
            for off in (
                SVector(-r, -r, -r), SVector(r, -r, -r), SVector(-r, r, -r), SVector(r, r, -r),
                SVector(-r, -r, r), SVector(r, -r, r), SVector(-r, r, r), SVector(r, r, r),
                SVector(r, 0.0, 0.0), SVector(-r, 0.0, 0.0), SVector(0.0, r, 0.0),
                SVector(0.0, -r, 0.0), SVector(0.0, 0.0, r), SVector(0.0, 0.0, -r),
            )
                x_try = old + off
                valid, cmin, cmax, sin2 = eval_candidate(mesh, v, x_try, elems, true)
                (valid && quality_not_worse(cmin, cmax, cmin_cur, cmax_cur)) || continue
                dmove = norm(x_try - old)
                if sin2 > best_sin2 || (sin2 == best_sin2 && dmove < best_d)
                    best_sin2 = sin2; best_x = x_try; best_d = dmove
                end
            end
        end
        if best_d > 0.0
            mesh.X[v] = best_x
            mesh.node_sdf[v] = eval_sdf(mesh, best_x)
            return (true, best_sin2)
        end
        return (false, sin2_cur)
    end
end

# ------------------------------------------------------------------------------
# Driver
# ------------------------------------------------------------------------------

"""
    relax_mesh!(mesh::BlockMesh, opts::RelaxOptions) -> mesh

Optional, gated vertex-relaxation post-pass (Etapa 10), run after the final
`update_connectivity!` so `mesh.INE` and the boundary triangulation are valid. Moves
vertices only -- never changes topology -- so the watertight trim structure,
determinism, and every node index survive. Surface vertices stay exactly on `phi = 0`;
interior vertices stay strictly inside; the local minimum quality never decreases, so
the global minimum dihedral is guaranteed not to drop. See [`RelaxOptions`](@ref).
"""
function relax_mesh!(mesh::BlockMesh, opts::RelaxOptions)
    isempty(mesh.IEN) && return mesh
    @info "Relaxing mesh (mode = :$(opts.mode), sizing = :$(opts.sizing), band = $(opts.band))..."

    infra = build_relax_infra(mesh, opts)
    if isempty(infra.movable)
        @info "  No movable vertices in the active band -- nothing to relax."
        return mesh
    end
    println("  Active band: $(length(infra.movable)) movable vertices " *
            "($(count(infra.is_surface[v] for v in infra.movable)) on the surface)")

    prev_quality = -Inf
    for sweep = 1:opts.max_sweeps
        moved = 0
        max_disp = 0.0
        sweep_quality = Inf
        for v in infra.movable
            elems = mesh.INE[v]
            if opts.mode === :uniform
                did_move, disp = relax_vertex_uniform!(mesh, v, infra, opts, elems)
                did_move && (moved += 1; max_disp = max(max_disp, disp))
            else
                did_move, qv = relax_vertex_quality!(mesh, v, infra, opts, elems)
                did_move && (moved += 1)
                sweep_quality = min(sweep_quality, qv)
            end
        end

        if opts.mode === :uniform
            @info "  sweep $sweep: moved $moved/$(length(infra.movable)), " *
                  "max displacement $(round(max_disp; digits = 5))"
            if max_disp < 1e-2 * mesh.grid_step
                @info "  Converged: displacement below tolerance."
                break
            end
        else
            @info "  sweep $sweep: moved $moved/$(length(infra.movable)), " *
                  "min quality (sin^2) $(round(sweep_quality; digits = 6))"
            if sweep_quality <= prev_quality
                @info "  Converged: quality no longer improving."
                break
            end
            prev_quality = sweep_quality
        end
    end

    # Bookkeeping: node_hash / node_map are coordinate-keyed connectivity-merge artifacts
    # (built by update_connectivity!). After moving nodes they are stale, and nothing in
    # the pipeline reads them past this point (export and the optional plane warp do not),
    # so empty them to avoid any later code trusting a stale coordinate -> index map.
    empty!(mesh.node_hash)
    empty!(mesh.node_map)
    return mesh
end
