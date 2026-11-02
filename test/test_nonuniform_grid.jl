# ==============================================================================
# Regression test: A15 tile nodes are keyed by their integer lattice coordinate
# ==============================================================================
# On a Float32 grid with a non-dyadic step (the smoothed-field files of rho2sdf.jl before
# its Float64 export, e.g. the gripper test data) the grid spacing varies by a few ulp.
# process_cell_A15! keyed the tile nodes by their position rounded to 1e-8 h, so a node
# shared by two cells was stored twice; warp! then moved the two copies apart and the
# final mesh had crack faces inside the solid. The cracks are closed surfaces, so
# check_watertight still reported 0 open edges (gripper: 34,513 such faces).
#
# Run stand-alone: julia --project=. test/test_nonuniform_grid.jl
# ==============================================================================
using Test, LinearAlgebra, StaticArrays
using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh

if !isdefined(@__MODULE__, :check_watertight)
    include(joinpath(@__DIR__, "helpers", "helpers.jl"))
end

# Number of node pairs closer than r (hash grid of cell size r).
function _close_pairs(X, r)
    key(p) = Tuple(floor.(Int, p ./ r))
    cells = Dict{NTuple{3,Int},Vector{Int}}()
    for (i, p) in enumerate(X)
        push!(get!(cells, key(p), Int[]), i)
    end
    n = 0
    for (i, p) in enumerate(X), dz in -1:1, dy in -1:1, dx in -1:1
        k = key(p)
        for j in get(cells, (k[1] + dx, k[2] + dy, k[3] + dz), Int[])
            j > i && norm(X[j] - p) < r && (n += 1)
        end
    end
    return n
end

@testset "A15 nodes keyed by lattice coordinate (non-uniform Float32 grid)" begin
    # Sphere on a 25^3 grid with step 0.3 starting at 10, in the data-file convention
    # (positive inside). The Float32 axis is built like the old rho2sdf create_smooth_grid.
    h = 0.3; n = 25; o = 10.0
    c = SVector(o + 3.613, o + 3.621, o + 3.607); r = 2.4
    phi(p) = r - norm(SVector{3,Float64}(p...) - c)
    ax = range(Float32(o), Float32(o + (n - 1) * h), length = n)
    grid32 = [Float32[x, y, z] for x in ax, y in ax, z in ax]
    grid64 = [o .+ h .* [i, j, k] for i = 0:n-1, j = 0:n-1, k = 0:n-1]
    sdf32 = [Float32(phi(p)) for p in grid32]
    sdf64 = [phi(p) for p in grid64]

    # The Float32 spacing deviates from h far beyond the old 1e-8 h node tolerance.
    @test maximum(abs(ax[i+1] - ax[i] - h) for i = 1:n-1) > 1e-6 * h

    # Fill: every lattice node exactly once, the same nodes as on the exact Float64 grid
    # (before the fix: 23,686 nodes with 3,054 duplicated pairs vs 20,632).
    m32 = BlockMesh(sdf32, grid32)
    m64 = BlockMesh(sdf64, grid64)
    redirect_stdout(devnull) do
        generate_mesh!(m32, "A15")
        generate_mesh!(m64, "A15")
    end
    @test length(m32.X) == length(m64.X)
    @test _close_pairs(m32.X, 1e-3 * h) == 0

    # Final mesh: no duplicated nodes and no boundary faces inside the solid (before the
    # fix: 844 boundary faces with the centroid deeper than 0.1 h, 6 non-manifold edges).
    mesh = mktempdir() do d
        redirect_stdout(devnull) do
            generate_tetrahedral_mesh(BlockMesh(sdf32, grid32), joinpath(d, "sphere32"))
        end
    end
    w = check_watertight(mesh)
    @test w.open_edges == 0
    @test w.nonmanifold_edges == 0
    @test _close_pairs(mesh.X, 1e-3 * h) == 0
    faces = Dict{NTuple{3,Int},Int}()
    for t in mesh.IEN, f in tet_faces(t)
        faces[f] = get(faces, f, 0) + 1
    end
    @test all(f -> phi((mesh.X[f[1]] + mesh.X[f[2]] + mesh.X[f[3]]) / 3) < 0.1 * h,
              [f for (f, k) in faces if k == 1])
end
