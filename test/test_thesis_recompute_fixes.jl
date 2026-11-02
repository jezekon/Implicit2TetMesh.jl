# ==============================================================================
# Regression tests for the session-5 fixes of the dissertation recompute (commit
# b1985a2, tag thesis-recompute-s5). Each failed before its fix.
#
#   C-F1  process_cell_A15! skipped a cell when its 8 grid corners and those of its 6
#         face neighbours were outside, although the A15 tile has nodes at h/4 inside
#         the cell: a plate thinner than h between two grid planes gave an empty mesh.
#   C-F2  inverse_hex_map reported inside = true for a Newton iterate that had not
#         converged (residual of several mm), so UnstructuredSDF located points in
#         elements that do not contain them.
#   C-F3  the node alignment projected a node onto its plane even when that inverted an
#         incident tet; remove_inverted_elements! then reoriented the tet (a fold).
#
# Run stand-alone: julia --project=. test/test_thesis_recompute_fixes.jl
# ==============================================================================
using Test, LinearAlgebra, StaticArrays
using Implicit2TetMesh
import Implicit2TetMesh.Fundamentals: SDFSource, eval_sdf, bbox, inverse_hex_map, hex8_map_and_jacobian, BlockMesh
import Implicit2TetMesh.GenerateMesh: create_INE!
import Implicit2TetMesh.Modification: warp_mesh_by_planes_sdf!

# Analytic field as an SDF source (phi < 0 inside).
struct _AnalyticSDF{F} <: SDFSource
    f::F
    bmin::SVector{3,Float64}
    bmax::SVector{3,Float64}
end
bbox(s::_AnalyticSDF) = (s.bmin, s.bmax)
eval_sdf(s::_AnalyticSDF, p::SVector{3,Float64}) = s.f(p)

@testset "Session-5 fixes (thesis recompute)" begin

    @testset "C-F1: plate thinner than h between two grid planes is meshed" begin
        # Plate |z - 0.5| <= 0.3 with a rounded rim (points within 0.3 of a flat disk of
        # radius 5). With dx = 1 and this bbox the grid planes are z = ..., 0, 1, ...: no grid
        # corner lies inside the plate, the A15 tile nodes at z = 0.25, 0.5, 0.75 do.
        c = SVector(0.11, 0.07, 0.5); r = 0.3; R = 5.0
        f = p -> begin
            q = p - c; z = q[3]; ρ = hypot(q[1], q[2])
            (ρ <= R ? abs(z) : sqrt((ρ - R)^2 + z^2)) - r
        end
        src = _AnalyticSDF(f, SVector(-8.0, -8.0, -8.0), SVector(8.0, 8.0, 8.0))
        opts = MeshGenerationOptions(cut_points = :bisection)
        mesh = mktempdir() do d
            cd(d) do
                redirect_stdout(devnull) do
                    generate_tetrahedral_mesh(src, "plate"; dx = 1.0, options = opts)
                end
            end
        end
        @test length(mesh.IEN) > 0
        vol = sum(abs(dot(mesh.X[t[2]] - mesh.X[t[1]], cross(mesh.X[t[3]] - mesh.X[t[1]], mesh.X[t[4]] - mesh.X[t[1]]))) / 6
                  for t in mesh.IEN; init = 0.0)
        vexact = 2π * R^2 * r + π^2 * R * r^2 + 4 / 3 * π * r^3
        @test abs(vol - vexact) / vexact < 0.1
        @test all(i -> f(mesh.X[i]) <= 1e-9, unique(i for t in mesh.IEN for i in t))
    end

    @testset "C-F2: inverse_hex_map does not report a non-converged point as inside" begin
        # Hex 11972 of the gripper SIMP mesh (VTK order) and a point about 3 mm outside it,
        # for which the old code returned inside = true with a residual of 6.2 mm.
        X = (SVector(15.036155700684001, -45.043384552002003, 105.0),
             SVector(16.129999160766999, -41.259998321532997, 105.0),
             SVector(12.050000190735, -42.25, 105.0),
             SVector(10.609999656676999, -45.169998168945, 105.0),
             SVector(15.036155700684001, -45.043384552002003, 110.0),
             SVector(16.129999160766999, -41.259998321532997, 110.0),
             SVector(12.050000190735, -42.25, 110.0),
             SVector(10.609999656676999, -45.169998168945, 110.0))
        q = SVector(10.067699999999999, -38.112099999999998, 106.71899999999999)
        _, inside = inverse_hex_map(X, q)
        @test !inside
        # A point inside the same hex is still found, with its natural coordinates.
        p, _ = hex8_map_and_jacobian(X, SVector(0.3, -0.4, 0.5))
        xi, inside2 = inverse_hex_map(X, p)
        @test inside2
        @test isapprox(xi, SVector(0.3, -0.4, 0.5); atol = 1e-9)
    end

    @testset "C-F3: the node alignment does not invert a tet" begin
        # One tet ABCD, D at z = 0.1 above the face ABC (z = 0). A small plane patch at
        # z = -0.15 under D (within the threshold 0.3 h of D, but A, B, C project outside the
        # patch): projecting D onto it would put D below ABC, i.e. invert the tet. A patch at
        # z = 0.2 is a harmless projection and must still be applied.
        src = _AnalyticSDF(p -> -1.0, SVector(0.0, 0.0, 0.0), SVector(1.0, 1.0, 1.0))
        function onetet()
            m = BlockMesh(src; dx = 1.0, padding = 0)
            m.X = [SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0), SVector(0.0, 1.0, 0.0), SVector(0.3, 0.3, 0.1)]
            m.IEN = [[1, 2, 3, 4]]
            m.node_sdf = zeros(4)
            create_INE!(m)
            return m
        end
        m = onetet()
        bad = [PlaneDefinition([0.0, 0.0, 1.0], [0.3, 0.3, -0.15], Circle(0.05))]
        redirect_stdout(devnull) do
            warp_mesh_by_planes_sdf!(m, bad, 0.3)
        end
        @test m.X[4] == SVector(0.3, 0.3, 0.1)          # projection rejected, node left
        @test m.IEN == [[1, 2, 3, 4]]                   # tet neither reoriented nor removed
        m = onetet()
        good = [PlaneDefinition([0.0, 0.0, 1.0], [0.3, 0.3, 0.2], Circle(0.05))]
        redirect_stdout(devnull) do
            warp_mesh_by_planes_sdf!(m, good, 0.3)
        end
        @test abs(m.X[4][3] - 0.2) < 1e-5               # harmless projection applied
        @test m.IEN == [[1, 2, 3, 4]]
    end
end
