# Cross-validation of Implicit2TetMesh against quartet (Doran, Chang & Bridson 2013,
# github.com/crawforddoran/quartet), the reference implementation of isosurface stuffing on
# the A15 lattice. Not part of Pkg.test(): it needs a C++ build of quartet. See
# doc/VERIFICATION.md for the method and the recorded results.
#
# Both programs get bit-identical input: the package's trilinear field sampled on the input
# grid refined 4x (where quartet's 4-cell tile and the package's 1-cell tile coincide),
# rounded to Float32 (quartet's precision). quartet runs through quartet_driver.cpp
# (unmodified quartet code; dumps the mesh after cut_lattice / warp_vertices / trim_spikes /
# remove_exterior_tets); the package runs the same steps on a StructuredSDF of the same values.
#
#   git clone https://github.com/crawforddoran/quartet.git <Q>
#   clang++ -O2 -std=c++11 -w -I<Q>/src -o test/quartet_crossval/quartet_driver \
#       test/quartet_crossval/quartet_driver.cpp <Q>/src/{predicates,geometry_queries,sdf,trimesh,tet_mesh,feature,read_obj,tet_quality,match_features,optimize_tet_mesh,make_signed_distance}.cpp
#   julia --project=. test/quartet_crossval/xval.jl beam <work dir>      # or: gripper
using Implicit2TetMesh, JLD2, StaticArrays, LinearAlgebra, Printf
import Implicit2TetMesh.Fundamentals: BlockMesh, StructuredSDF, eval_sdf
import Implicit2TetMesh.GenerateMesh: generate_mesh!, warp!, update_connectivity!,
    slice_ambiguous_tetrahedra!, remove_inverted_elements!, is_positively_oriented,
    compute_dihedral_angle_range
import Implicit2TetMesh.Modification: remove_isolated_components!
const STATE = normpath(joinpath(@__DIR__, "..", ".."))   # repository root
const CASE, WORK = ARGS[1], ARGS[2]
mkpath(WORK)
files = CASE == "beam" ?
    ("beam/Z_beam_HEX8_FineGrid_B-1.0_smooth-1.jld2", "beam/Z_beam_HEX8_FineSDF_B-1.0_smooth-1.jld2") :
    ("gripper/Z_robot_gripper_HEX8_FineGrid_B-2.0085_smooth-1.jld2", "gripper/Z_robot_gripper_HEX8_FineSDF_B-2.0085_smooth-1.jld2")
@load joinpath(STATE, "data", files[1]) fine_grid
@load joinpath(STATE, "data", files[2]) fine_sdf
coarse = BlockMesh(fine_sdf, fine_grid)              # standard structured path (trilinear source)
o = coarse.grid[1, 1, 1]; h = coarse.grid[2, 1, 1][1] - o[1]
n = size(coarse.grid); N = 4 .* (n .- 1) .+ 1
@printf("%s: input grid %s, origin %s, h = %.17g; refined grid %s\n", CASE, n, o, h, N)

# ---- identical Float32 input for both
hf = h / 4
fgrid = [SVector{3,Float64}(o[1] + hf * (a - 1), o[2] + hf * (b - 1), o[3] + hf * (c - 1))
         for a in 1:N[1], b in 1:N[2], c in 1:N[3]]
phi32 = Float32[eval_sdf(coarse.sdf_source, p) for p in fgrid]
open(joinpath(WORK, "sdf.bin"), "w") do io
    write(io, Int32.(collect(N))); write(io, Float32.([o[1], o[2], o[3], hf])); write(io, phi32)
end
# float32 grid coordinates as quartet computes them: origin + dx*Vec3f(vertex)
@assert all(Float32(o[d]) == o[d] for d in 1:3) && Float32(hf) == hf "origin/spacing not exact in Float32"
run(pipeline(`$(joinpath(@__DIR__, "quartet_driver")) $(joinpath(WORK, "sdf.bin")) $(joinpath(WORK, "q"))`; stdout = stdout))

function readq(step)
    open(joinpath(WORK, "q_$step.bin")) do io
        nv, nt = read(io, Int32), read(io, Int32)
        X = reshape(read!(io, Vector{Float32}(undef, 3nv)), 3, nv)
        vphi = read!(io, Vector{Float32}(undef, nv))
        T = reshape(read!(io, Vector{Int32}(undef, 4nt)), 4, nt) .+ 1
        (X = Float64.(X), vphi = vphi, T = T)
    end
end

# ---- package on the same Float32 values
src = StructuredSDF(fgrid, Float64.(phi32))
mesh = BlockMesh(src; dx = h, padding = 0)
@assert mesh.grid[1, 1, 1] == o && size(mesh.grid) == n "lattice not aligned with quartet's tiles"
quiet(f) = redirect_stdout(f, devnull)

# vertex matching by position (tolerance tol·h): quartet computes positions in Float32
# (origin + dx*i), so they differ from the package's Float64 ones by ~1e-6 h for a
# non-dyadic spacing; distinct vertices are >= 0.05 h apart, so 1e-4 h is unambiguous
function matcher(Xq, tol)
    key(p) = Tuple(round.(Int, p ./ (tol * h)))
    d = Dict{NTuple{3,Int},Vector{Int}}()
    for i in axes(Xq, 2); push!(get!(d, key(Xq[:, i]), Int[]), i); end
    function find(p)
        k = key(p); best = 0; bd = Inf
        for dz in -1:1, dy in -1:1, dx in -1:1
            for i in get(d, (k[1] + dx, k[2] + dy, k[3] + dz), Int[])
                dd = norm(Xq[:, i] - p); dd < bd && (bd = dd; best = i)
            end
        end
        return bd <= tol * h ? best : 0
    end
    return find
end
tetset(T) = Set(Tuple(sort(T[:, t])) for t in axes(T, 2))
mapped(ien, mp) = [Tuple(sort([mp[v] for v in t])) for t in ien]

report = String[]
say(s) = (println(s); push!(report, s))

# 1. lattice
quiet(() -> generate_mesh!(mesh, "A15"))
q1 = readq("1_lattice")
find1 = matcher(q1.X, 1e-4)
mp = [find1(p) for p in mesh.X]
all(>(0), mp) || error("$(count(==(0), mp)) package nodes have no quartet counterpart")
say(@sprintf("1 lattice : package %d nodes %d tets | quartet %d vertices %d tets | unmatched nodes %d",
             length(mesh.X), length(mesh.IEN), size(q1.X, 2), size(q1.T, 2), count(==(0), mp)))
s_pkg = Set(mapped(mesh.IEN, mp)); s_q = tetset(q1.T)
dphi = [abs(mesh.node_sdf[i] - q1.vphi[mp[i]]) for i in eachindex(mesh.X)]
say(@sprintf("            tet sets equal: %s (only package %d, only quartet %d); node phi: %d of %d bitwise equal, max |diff| = %.1e h, same sign everywhere: %s",
             s_pkg == s_q, length(setdiff(s_pkg, s_q)), length(setdiff(s_q, s_pkg)), count(iszero, dphi), length(dphi),
             maximum(dphi) / h, all(i -> sign(mesh.node_sdf[i]) == sign(q1.vphi[mp[i]]), eachindex(mesh.X))))

# 2. warp (same vertex identities as step 1)
quiet(() -> warp!(mesh, "A15"))
q2 = readq("2_warp")
w_pkg = Set(i for i in eachindex(mesh.X) if mesh.node_sdf[i] == 0.0)
w_q = Set(i for i in eachindex(mesh.X) if q2.vphi[mp[i]] == 0)
dx = [norm(mesh.X[i] - q2.X[:, mp[i]]) for i in eachindex(mesh.X)]
# A vertex may have two candidate cut points whose metrics alpha*|e|^2 differ by less than
# quartet's Float32 resolution; quartet may then pick the other one (a near tie, not a rule
# difference). Such vertices are reported and their tets left out of the step-3 comparison.
moved = [i for i in eachindex(mesh.X) if dx[i] > 1e-4 * h]
say(@sprintf("2 warp    : warped package %d, quartet %d, same vertices: %s; max |x_pkg - x_quartet| = %.2e h over the %d others",
             length(w_pkg), length(w_q), w_pkg == w_q, maximum(dx[i] for i in eachindex(dx) if !(i in moved); init = 0.0) / h,
             length(mesh.X) - length(moved)))
say(@sprintf("            warped to a different cut point (Float32 near tie in quartet): %d vertices", length(moved)))
moved_pos = Set(mesh.X[i] for i in moved)
moved_q = Set(mp[i] for i in moved)

# 3. slice = trim_spikes + quadruple-zero resolution; quartet: trim_spikes, remove_exterior_tets
nX0 = length(mesh.X)
quiet(() -> update_connectivity!(mesh; build_ine = false))
say(@sprintf("  merge   : update_connectivity! after warp: %d -> %d nodes", nX0, length(mesh.X)))
quiet(() -> slice_ambiguous_tetrahedra!(mesh, "A15"))
q3 = readq("3_trim"); q4 = readq("4_exterior")
find3 = matcher(q3.X, 1e-4)
mp3 = [find3(p) for p in mesh.X]
say(@sprintf("3 slice   : package %d tets | quartet trim_spikes %d, after remove_exterior_tets %d | unmatched package nodes %d",
             length(mesh.IEN), size(q3.T, 2), size(q4.T, 2), count(==(0), mp3)))
# leave out the tets around the near-tie vertices: on the package side every tet with a
# moved vertex or an unmatched node (a cut point on an edge of a moved vertex), on the
# quartet side every tet with its moved vertex or a cut vertex the package does not have
pkg_bad(t) = any(v -> mp3[v] == 0 || mesh.X[v] in moved_pos, t)
used_q = Set(v for v in mp3 if v > 0)
q_bad(t) = any(v -> v in moved_q || !(v in used_q), t)
n_ex_p = count(pkg_bad, mesh.IEN); n_ex_q = count(t -> q_bad(q4.T[:, t]), axes(q4.T, 2))
n_ex_p + n_ex_q > 0 && say(@sprintf("            left out around the near-tie vertices: package %d tets, quartet %d tets", n_ex_p, n_ex_q))
s_pkg = Set(mapped([t for t in mesh.IEN if !pkg_bad(t)], mp3))
s_q4 = Set(Tuple(sort(q4.T[:, t])) for t in axes(q4.T, 2) if !q_bad(q4.T[:, t]))
s_q3 = tetset(q3.T)
qz(t) = all(v -> q3.vphi[v] == 0, t)           # quartet quadruple-zero tet
only_pkg = setdiff(s_pkg, s_q4); only_q = setdiff(s_q4, s_pkg)
say(@sprintf("            non-quadruple-zero tets equal: %s",
             Set(t for t in s_pkg if !qz(t)) == Set(t for t in s_q4 if !qz(t))))
say(@sprintf("            quadruple-zero: package %d, quartet %d | only package %d, only quartet %d",
             count(qz, s_pkg), count(qz, s_q4), length(only_pkg), length(only_q)))
# classify the differences by the documented rules (dihedral filter [10, 140] deg, pocket rule)
inv_mp3 = Dict(v => i for (i, v) in enumerate(mp3) if v > 0)
function classify(t)
    tp = [inv_mp3[v] for v in t]
    c = sum(q3.X[:, v] for v in t) / 4
    centroid_out = eval_sdf(src, SVector{3,Float64}(c)) > 0
    tet = collect(tp)
    is_positively_oriented(mesh, tet) || (tet = [tet[1], tet[2], tet[4], tet[3]])
    lo, hi = compute_dihedral_angle_range(mesh, tet)
    return (centroid_out = centroid_out, bad_shape = lo < 10 || hi > 140, in_trim = t in s_q3)
end
cq = [classify(t) for t in only_q]; cp = [classify(t) for t in only_pkg]
say(@sprintf("            only quartet: %d, of them rejected by the [10,140] deg shape filter: %d",
             length(cq), count(c -> c.bad_shape, cq)))
say(@sprintf("            only package: %d, of them centroid outside (kept by the pocket rule): %d, in quartet's trim output: %d",
             length(cp), count(c -> c.centroid_out, cp), count(c -> c.in_trim, cp)))

# 4. rest of the package pipeline (quartet: compactMesh only)
quiet() do
    update_connectivity!(mesh; build_ine = false); remove_inverted_elements!(mesh)
    remove_isolated_components!(mesh, keep_largest = true); update_connectivity!(mesh)
end
say(@sprintf("4 final   : package %d nodes %d tets (after remove_inverted/remove_isolated) | quartet %d tets",
             length(mesh.X), length(mesh.IEN), size(q4.T, 2)))
write(joinpath(WORK, "report.txt"), join(report, "\n") * "\n")
