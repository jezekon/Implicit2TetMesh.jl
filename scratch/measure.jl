# Regression harness for the staged Implicit2TetMesh.jl work (see doc/PROMPTS.md).
# Run from the repo root, e.g.:  julia --project=. scratch/measure.jl beam
#
# The dihedral statistics are purely geometric (sign-agnostic), so they are valid
# in either SDF convention and can be compared before/after Etapa 1.
using Implicit2TetMesh, Implicit2TetMesh.Fundamentals, Implicit2TetMesh.GenerateMesh, Implicit2TetMesh.Modification
using JLD2, LinearAlgebra

nm = isempty(ARGS) ? "beam" : ARGS[1]
f = Dict(
    "beam" => ("Z_beam_HEX8_FineGrid_B-1.0_smooth-1.jld2", "Z_beam_HEX8_FineSDF_B-1.0_smooth-1.jld2"),
    "gripper" => ("Z_robot_gripper_HEX8_FineGrid_B-2.0085_smooth-1.jld2", "Z_robot_gripper_HEX8_FineSDF_B-2.0085_smooth-1.jld2"),
)
gf, sf = f[nm]
@load joinpath("data", nm, gf) fine_grid
@load joinpath("data", nm, sf) fine_sdf

# Dihedral-angle statistics + inverted-element count for a mesh.
function dih(mesh)
    mn = 180.0; mx = 0.0; l5 = 0; l10 = 0; g140 = 0; inv = 0
    for t in mesh.IEN
        v = [mesh.X[t[k]] for k = 1:4]
        dot(v[2] - v[1], cross(v[3] - v[1], v[4] - v[1])) <= 0 && (inv += 1)
        F = [(1, 2, 3, 4), (1, 2, 4, 3), (1, 3, 4, 2), (2, 3, 4, 1)]
        N = []; ok = true
        for (a, b, c, o) in F
            cr = cross(v[b] - v[a], v[c] - v[a]); n = norm(cr)
            n < 1e-14 && (ok = false; break)
            u = cr / n; ct = (v[a] + v[b] + v[c]) / 3
            dot(u, v[o] - ct) > 0 && (u = -u)
            push!(N, u)
        end
        ok || continue
        for i = 1:3, j = (i + 1):4
            a = acos(-clamp(dot(N[i], N[j]), -1, 1)) * 180 / pi
            mn = min(mn, a); mx = max(mx, a); a < 5 && (l5 += 1); a < 10 && (l10 += 1); a > 140 && (g140 += 1)
        end
    end
    (mn, mx, l5, l10, g140, inv)
end

m = BlockMesh(fine_sdf, fine_grid)
generate_mesh!(m, "A15"); warp!(m, "A15"); update_connectivity!(m)
slice_ambiguous_tetrahedra!(m, "A15"); update_connectivity!(m)
remove_inverted_elements!(m); remove_isolated_components!(m, keep_largest = true); update_connectivity!(m)
s = dih(m)
println("$nm: tets=$(length(m.IEN)) min=$(round(s[1];digits=3)) max=$(round(s[2];digits=3)) <5=$(s[3]) <10=$(s[4]) >140=$(s[5]) inv=$(s[6])")

# Optional Etapa-10 relaxation: pass a mode as the 2nd arg, e.g.
#   julia --project=. scratch/measure.jl beam uniform   (or: quality)
if length(ARGS) >= 2
    relax_mesh!(m, RelaxOptions(mode = Symbol(ARGS[2])))
    r = dih(m)
    println("$nm (relax=:$(ARGS[2])): tets=$(length(m.IEN)) min=$(round(r[1];digits=3)) max=$(round(r[2];digits=3)) <5=$(r[3]) <10=$(r[4]) >140=$(r[5]) inv=$(r[6])")
end
