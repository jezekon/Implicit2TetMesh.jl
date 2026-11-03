# ==============================================================================
# Domain clip of rho2sdf.jl: DomainDistance, ClippedRBFSDF and load_rbf_source
# ==============================================================================
# rho2sdf.jl (option domain_clip, on by default) clips its smoothed field by the design
# domain D, min(f - clip_shift, d_D), and writes clip_shift and the hex mesh of D into the
# JLD2. load_rbf_source then builds a ClippedRBFSDF, max(eval_sdf(rbf), -d_D) here.
#
#   (1) d_D against the exact signed distance of a box, on a block whose interior nodes are
#       jittered: the boundary faces stay planar, so d_D must still be exact, and the sign
#       test runs on distorted elements.
#   (2) load_rbf_source on a synthetic file: the source follows the key domain_clip,
#       domain_clip = false gives the unclipped RBFSDF with `shift`, and a file without the
#       clip loads as before.
#   (3) Mesh of a bar that runs through the whole RBF grid, clipped by D: the mesh ends at
#       the faces of D (the unclipped one does not) and passes the invariants.
#
# The values against rho2sdf.jl's own d_D and clipped field (cantilever, gripper) were
# checked outside the test suite: identical at 200,000 points each.
#
# Run stand-alone: julia --project=. test/test_clipped_rbf.jl
# ==============================================================================
using Test, LinearAlgebra, StaticArrays, Random, JLD2
using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals

if !isdefined(@__MODULE__, :check_watertight)
    include(joinpath(@__DIR__, "helpers", "helpers.jl"))
end

# Design domain: the block [-1, 1] × [-0.5, 0.5]² of 8 × 4 × 4 HEX8 elements (VTK corner
# order), as rho2sdf.jl exports it: nodes 3 × n, elements 8 × m. Interior nodes are moved
# by up to `jitter` of the element size; the boundary nodes stay, so D is the box.
function clip_test_block(; jitter = 0.0, seed = 3)
    rng = MersenneTwister(seed)
    nx, ny, nz = 8, 4, 4
    node(i, j, k) = 1 + i + (nx + 1) * (j + (ny + 1) * k)
    nodes = zeros(3, (nx + 1) * (ny + 1) * (nz + 1))
    for k = 0:nz, j = 0:ny, i = 0:nx
        p = [-1.0 + 0.25 * i, -0.5 + 0.25 * j, -0.5 + 0.25 * k]
        if 0 < i < nx && 0 < j < ny && 0 < k < nz
            p .+= jitter * 0.25 .* (2 .* rand(rng, 3) .- 1)
        end
        nodes[:, node(i, j, k)] = p
    end
    elements = zeros(Int, 8, nx * ny * nz)
    e = 0
    for k = 0:nz-1, j = 0:ny-1, i = 0:nx-1
        e += 1
        elements[:, e] = [node(i, j, k), node(i + 1, j, k), node(i + 1, j + 1, k), node(i, j + 1, k),
                          node(i, j, k + 1), node(i + 1, j, k + 1), node(i + 1, j + 1, k + 1), node(i, j + 1, k + 1)]
    end
    return nodes, elements
end

# Exact signed distance of the box [low, high], positive inside
function box_distance(p, low, high)
    gap = norm(max.(low .- p, 0.0, p .- high))
    return gap > 0 ? -gap : minimum(min.(p .- low, high .- p))
end

# A JLD2 as rho2sdf.jl writes it. RBF: a bar along x, |y|, |z| < 0.25, on the grid
# [-1.5, 1.5] × [-1, 1]², h = 0.125, which reaches beyond D at both ends. Weights of the
# normalized RBF approximation, w = g / S, S the sum of the kernel over the grid.
function write_clip_test_file(file; clip::Bool)
    h = 0.125
    origin = [-1.5, -1.0, -1.0]
    dims = (25, 17, 17)
    S = sum(exp(-(i^2 + j^2 + k^2)) for i = -2:2, j = -2:2, k = -2:2 if exp(-(i^2 + j^2 + k^2)) > 1e-3)
    weights = [(0.25 - max(abs(origin[2] + h * (j - 1)), abs(origin[3] + h * (k - 1)))) / S
               for i = 1:dims[1], j = 1:dims[2], k = 1:dims[3]]
    fields = (weights = weights, origin = origin, h = h, sigma = h, threshold = 1e-3, shift = 0.01,
              inside = "positive", smoothing = "approximation", target_volume = 0.0, volume = 0.0)
    if clip
        nodes, elements = clip_test_block()
        fields = merge(fields, (domain_clip = true, clip_shift = -0.02, clip_volume = 0.0,
                                volume_outside = 0.0, domain_nodes = nodes, domain_elements = elements))
    end
    jldsave(file; fields...)
    return file
end

@testset "Domain clip (rho2sdf.jl)" begin

    @testset "DomainDistance: exact on a box with distorted elements" begin
        cap = 0.4
        low, high = [-1.0, -0.5, -0.5], [1.0, 0.5, 0.5]
        for jitter in (0.0, 0.3)
            d_D = DomainDistance(clip_test_block(; jitter)..., cap)
            @test length(d_D.triangles) == 4 * 2 * (8 * 4 + 8 * 4 + 4 * 4)
            largest_error = 0.0
            for x in range(-1.3, 1.3, 23), y in range(-0.8, 0.8, 15),
                z in (-0.7, -0.5, -0.31, 0.0, 0.17, 0.5, 0.66)
                p = SVector(x, y, z)
                largest_error = max(largest_error, abs(d_D(p) - clamp(box_distance(p, low, high), -cap, cap)))
            end
            @test largest_error <= 1e-12
        end
        @test_throws ErrorException DomainDistance(zeros(2, 4), ones(Int, 8, 1), 1.0)
        @test_throws ErrorException DomainDistance(zeros(3, 4), fill(5, 8, 1), 1.0)
    end

    @testset "load_rbf_source follows the file" begin
        dir = mktempdir()
        clipped_file = write_clip_test_file(joinpath(dir, "clipped.jld2"); clip = true)
        plain_file = write_clip_test_file(joinpath(dir, "plain.jld2"); clip = false)

        src = load_rbf_source(clipped_file)
        @test src isa ClippedRBFSDF
        @test src.rbf.shift == -0.02                     # clip_shift, not shift
        @test src.h == 0.125
        @test bbox(src) == bbox(src.rbf)

        rng = MersenneTwister(5)
        low, high = bbox(src)
        matches = true
        for _ = 1:2000
            p = low .+ rand(rng, SVector{3,Float64}) .* (high .- low)
            matches &= eval_sdf(src, p) == max(eval_sdf(src.rbf, p), -src.domain(p))
        end
        @test matches
        @test eval_sdf(src, SVector(5.0, 0.0, 0.0)) == eval_sdf(src.rbf, SVector(5.0, 0.0, 0.0))

        unclipped = load_rbf_source(clipped_file; domain_clip = false)
        @test unclipped isa RBFSDF
        @test unclipped.shift == 0.01

        @test load_rbf_source(plain_file) isa RBFSDF
        @test_throws ErrorException load_rbf_source(plain_file; domain_clip = true)
    end

    @testset "Mesh of the clipped bar ends at the domain" begin
        dir = mktempdir()
        file = write_clip_test_file(joinpath(dir, "clipped.jld2"); clip = true)
        options = MeshGenerationOptions(cut_points = :bisection)
        meshes = Dict{Bool,Any}()
        for clip in (true, false)
            src = load_rbf_source(file; domain_clip = clip)
            meshes[clip] = cd(() -> generate_tetrahedral_mesh(src, "bar"; dx = src.h, options = options), dir)
        end
        m = meshes[true]
        x = [p[1] for p in m.X]

        @test check_watertight(m).open_edges == 0
        @test count_inverted_exact(m) == 0
        @test count_components(m) == 1
        @test min_signed_volume(m) > 0

        # Within the domain up to the bisection tolerance, with flat ends on x = ±1
        @test maximum(abs.(x)) <= 1 + 1e-9
        @test count(abs.(x) .>= 1 - 1e-9) > 0
        # The unclipped bar runs on beyond the domain
        @test maximum(abs(p[1]) for p in meshes[false].X) > 1.2
    end
end
