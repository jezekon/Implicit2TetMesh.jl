# ==============================================================================
# Etapa 8 -- pluggable SDF sources (structured trilinear + unstructured HEX8 FE)
# ==============================================================================
# Three layers, from cheap to full-pipeline:
#
#   (1) UNIT    -- the source machinery in isolation: HEX8 shape functions reproduce
#       trilinear interpolation on a cube, the inverse isoparametric (Newton) map
#       round-trips on a distorted hex, the auto-detector recognises a uniform grid
#       and rejects a jittered one, the SIMP / level-set adapters honour their
#       conventions, and eval_sdf on a structured source equals eval_sdf on the SAME
#       field re-expressed as HEX8 elements (plus the outside-domain rule).
#
#   (2) ROUND-TRIP -- re-express the real beam grid as a conforming HEX8 mesh, force
#       the unstructured code path, and require the SAME output as the structured beam
#       run within floating-point tolerance. quartet can only consume structured grids,
#       so there is NO oracle for unstructured inputs: the structured run IS the oracle
#       here. This validates point location + inverse map + shape functions on real data.
#
#   (3) ANALYTIC -- a graded + jittered HEX8 block carrying an analytic sphere field,
#       meshed through the FE path. Judged (no oracle) by watertightness, zero inverted
#       tets, a single component, a dihedral quality floor, geometric fidelity to the
#       TRUE sphere, and a frozen node/tet baseline (re-freezable).

using Test
using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh
using StaticArrays, LinearAlgebra, Random

const F = Implicit2TetMesh.Fundamentals     # reach the internal (unexported) helpers

# Idempotent include guard so this file also runs stand-alone.
if !isdefined(@__MODULE__, :check_watertight)
    include(joinpath(@__DIR__, "helpers", "helpers.jl"))
end
if !isdefined(@__MODULE__, :BEAM_NODES)
    include(joinpath(@__DIR__, "helpers", "baselines.jl"))
end

# ------------------------------------------------------------------------------
# Local mesh-construction helpers (specific to these tests)
# ------------------------------------------------------------------------------

# VTK HEX8 corner offsets for grid cell (i,j,k); the SAME order get_cell_sdf_values
# and HEX8_NATURAL use, so a grid re-expressed this way interpolates consistently.
const HEX_OFFSETS = (
    (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
    (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
)

"""
    grid_to_hex(grid, vals) -> (nodes, hexes, phi)

Re-express a structured grid (node coordinates + nodal field) as a conforming HEX8
mesh: one hex per grid cell, in VTK corner order. Used to push a known structured
field through the unstructured code path.
"""
function grid_to_hex(grid::Array{SVector{3,Float64},3}, vals::Array{Float64,3})
    nx, ny, nz = size(grid)
    lin = LinearIndices((nx, ny, nz))
    nodes = Vector{SVector{3,Float64}}(undef, nx * ny * nz)
    phi = Vector{Float64}(undef, nx * ny * nz)
    for k = 1:nz, j = 1:ny, i = 1:nx
        nodes[lin[i, j, k]] = grid[i, j, k]
        phi[lin[i, j, k]] = vals[i, j, k]
    end
    hexes = NTuple{8,Int}[]
    for k = 1:nz-1, j = 1:ny-1, i = 1:nx-1
        push!(hexes, ntuple(a -> (o = HEX_OFFSETS[a]; lin[i+o[1], j+o[2], k+o[3]]), 8))
    end
    return nodes, hexes, phi
end

"""
    sphere_in_distorted_block(; N, L, center, radius, jitter, seed) -> (nodes, hexes, phi)

A genuinely unstructured input: a conforming HEX8 block on `[0,L]^3` with GRADED 1-D
levels (non-uniform spacing defeats the structured auto-detector) and interior-node
JITTER (makes the elements non-affine, so the Newton inverse map is exercised).
Boundary nodes are left in place, so the domain box stays exactly `[0,L]^3`. The nodal
field is the analytic sphere SDF `|x-center|-radius` (already phi<0=inside).
"""
function sphere_in_distorted_block(; N = 12, L = 4.0, center = SVector(2.0, 2.0, 2.0),
                                   radius = 1.2, jitter = 0.15, seed = 7)
    rng = MersenneTwister(seed)
    levels() = begin
        w = [1.0 + 0.6 * sin(3.0 * i / N) for i = 0:N-1]   # positive, non-uniform spacings
        x = zeros(N + 1)
        for i = 1:N
            x[i+1] = x[i] + w[i]
        end
        x .*= L / x[end]                                    # rescale to span [0,L]
        x
    end
    xs, ys, zs = levels(), levels(), levels()
    np = N + 1
    lin = LinearIndices((np, np, np))
    nodes = Vector{SVector{3,Float64}}(undef, np^3)
    for k = 1:np, j = 1:np, i = 1:np
        p = SVector(xs[i], ys[j], zs[k])
        if (1 < i < np) && (1 < j < np) && (1 < k < np)
            h = min(xs[i+1] - xs[i], xs[i] - xs[i-1],
                    ys[j+1] - ys[j], ys[j] - ys[j-1],
                    zs[k+1] - zs[k], zs[k] - zs[k-1])
            p = p + jitter * h * SVector{3,Float64}((2 .* rand(rng, 3) .- 1)...)
        end
        nodes[lin[i, j, k]] = p
    end
    hexes = NTuple{8,Int}[]
    for k = 1:np-1, j = 1:np-1, i = 1:np-1
        push!(hexes, ntuple(a -> (o = HEX_OFFSETS[a]; lin[i+o[1], j+o[2], k+o[3]]), 8))
    end
    phi = [norm(p - center) - radius for p in nodes]
    return nodes, hexes, phi
end

# Indices of boundary vertices (vertices of faces incident to exactly one tet).
function boundary_vertex_indices(mesh::BlockMesh)
    face_count = Dict{NTuple{3,Int},Int}()
    for tet in mesh.IEN, f in tet_faces(tet)
        face_count[f] = get(face_count, f, 0) + 1
    end
    bset = Set{Int}()
    for (f, c) in face_count
        c == 1 && (push!(bset, f[1]); push!(bset, f[2]); push!(bset, f[3]))
    end
    return collect(bset)
end

@testset "SDF sources (Etapa 8)" begin

    # --------------------------------------------------------------------------
    @testset "Unit: shape functions, inverse map, adapters" begin
        Random.seed!(11)

        # HEX8 shape functions on a unit cube reproduce trilinear interpolation.
        Xc = ntuple(a -> SVector{3,Float64}(HEX_OFFSETS[a]...), 8)
        cvals = SVector{8,Float64}(randn(8))
        manual_tri(c, x, y, z) = begin
            c00 = c[1] * (1 - x) + c[2] * x;  c10 = c[4] * (1 - x) + c[3] * x
            c01 = c[5] * (1 - x) + c[6] * x;  c11 = c[8] * (1 - x) + c[7] * x
            (c00 * (1 - y) + c10 * y) * (1 - z) + (c01 * (1 - y) + c11 * y) * z
        end
        max_tri = 0.0
        for _ = 1:1000
            x, y, z = rand(), rand(), rand()
            xi, inside = F.inverse_hex_map(Xc, SVector(x, y, z))
            @test inside
            max_tri = max(max_tri, abs(dot(F.hex8_shape(xi), cvals) - manual_tri(cvals, x, y, z)))
        end
        @test max_tri < 1e-12

        # Inverse isoparametric map round-trips on a distorted (non-affine) hex.
        Xd = ntuple(a -> Xc[a] + 0.15 * SVector{3,Float64}(randn(3)...), 8)
        @test det(F.hex8_map_and_jacobian(Xd, zero(SVector{3,Float64}))[2]) > 0
        max_inv = 0.0
        for _ = 1:1000
            xi_true = SVector{3,Float64}((2 .* rand(3) .- 1)...)
            p, _ = F.hex8_map_and_jacobian(Xd, xi_true)
            xi_rec, inside = F.inverse_hex_map(Xd, p)
            @test inside
            max_inv = max(max_inv, maximum(abs.(xi_rec - xi_true)))
        end
        @test max_inv < 1e-9

        # Auto-detector: uniform grid -> StructuredSDF (fast path); jittered -> HEX8.
        nx, ny, nz = 5, 4, 3
        grid = Array{SVector{3,Float64},3}(undef, nx, ny, nz)
        vals = Array{Float64,3}(undef, nx, ny, nz)
        for k = 1:nz, j = 1:ny, i = 1:nx
            grid[i, j, k] = SVector{3,Float64}(2.0 * (i - 1), 1.5 * (j - 1), 0.7 * (k - 1))
            vals[i, j, k] = sin(i) + cos(j) + k
        end
        nodes, hexes, phi = grid_to_hex(grid, vals)
        @test build_sdf_source(nodes, hexes, phi) isa StructuredSDF
        @test build_sdf_source(nodes, hexes, phi; force = :unstructured) isa UnstructuredSDF
        @test_throws ErrorException build_sdf_source(nodes, hexes, phi; force = :bogus)
        jit = copy(nodes)
        lin = LinearIndices((nx, ny, nz))
        for k = 2:nz-1, j = 2:ny-1, i = 2:nx-1
            jit[lin[i, j, k]] += SVector{3,Float64}(0.1, -0.1, 0.05)
        end
        @test build_sdf_source(jit, hexes, phi) isa UnstructuredSDF
        @test_throws ErrorException build_sdf_source(jit, hexes, phi; force = :structured)

        # Structured vs unstructured eval_sdf agree on a grid-derived field (the key
        # invariant behind the round trip), and bbox + outside-domain rule are correct.
        gc = SVector(4.0, 3.0, 2.5)
        for k = 1:nz, j = 1:ny, i = 1:nx
            vals[i, j, k] = norm(grid[i, j, k] - gc) - 2.0
        end
        nodes, hexes, phi = grid_to_hex(grid, vals)
        s_struct = StructuredSDF(grid, vals)
        s_unstr = build_sdf_source(nodes, hexes, phi; force = :unstructured)
        (bmin, bmax) = bbox(s_unstr)
        @test bbox(s_struct) == (bmin, bmax)
        max_eval = 0.0
        for _ = 1:5000
            p = SVector{3,Float64}(
                bmin[1] + rand() * (bmax[1] - bmin[1]),
                bmin[2] + rand() * (bmax[2] - bmin[2]),
                bmin[3] + rand() * (bmax[3] - bmin[3]))
            max_eval = max(max_eval, abs(eval_sdf(s_struct, p) - eval_sdf(s_unstr, p)))
        end
        @test max_eval < 1e-9
        off = SVector(1.0, 2.0, 0.5)
        @test eval_sdf(s_unstr, bmax + off) ≈ norm(off)     # outside -> distance to box

        # SIMP adapter: element densities -> nodal phi = iso_level - rho (phi<0=solid).
        dens_all_solid = ones(length(hexes))
        s_simp = simp_to_sdf_source(nodes, hexes, dens_all_solid; iso_level = 0.5)
        cpt = (bmin + bmax) / 2
        @test eval_sdf(s_simp, cpt) ≈ -0.5                  # rho=1 everywhere -> 0.5-1
        s_void = simp_to_sdf_source(nodes, hexes, zeros(length(hexes)); iso_level = 0.5)
        @test eval_sdf(s_void, cpt) ≈ 0.5                   # rho=0 -> 0.5-0
        @test_throws ErrorException simp_to_sdf_source(nodes, hexes, ones(length(hexes) + 1))

        # Level-set adapter: sign convention maps the input to phi<0=inside.
        neg = phi                                            # already negative inside
        pos = -phi                                           # positive inside
        s_neg = levelset_to_sdf_source(nodes, hexes, neg; inside = :negative)
        s_pos = levelset_to_sdf_source(nodes, hexes, pos; inside = :positive)
        for _ = 1:200
            p = SVector{3,Float64}(
                bmin[1] + rand() * (bmax[1] - bmin[1]),
                bmin[2] + rand() * (bmax[2] - bmin[2]),
                bmin[3] + rand() * (bmax[3] - bmin[3]))
            @test eval_sdf(s_neg, p) ≈ eval_sdf(s_pos, p)
        end
        @test_throws ErrorException levelset_to_sdf_source(nodes, hexes, neg; inside = :sideways)
    end

    # --------------------------------------------------------------------------
    # The structured beam run is the oracle: the same field re-expressed as HEX8 and
    # forced through the FE path must reproduce it within floating-point tolerance.
    @testset "Round-trip: beam as unstructured HEX8" begin
        outdir = ensure_output_dir()
        fine_grid, fine_sdf = load_beam()

        # Structured reference (the oracle).
        mesh_s = BlockMesh(fine_sdf, fine_grid)

        # Re-express its lattice field (already negated to phi<0=inside) as HEX8 and
        # force the unstructured path.
        nodes, hexes, phi = grid_to_hex(mesh_s.grid, mesh_s.SDF)
        src = build_sdf_source(nodes, hexes, phi; force = :unstructured)
        @test src isa UnstructuredSDF
        @test bbox(src) == (mesh_s.grid[1, 1, 1], mesh_s.grid[end, end, end])

        # padding=0 reproduces the structured lattice exactly (integer * dx), and the
        # corner sampling reproduces the cached field bit-for-bit.
        mesh_u = BlockMesh(src; dx = mesh_s.grid_step, padding = 0)
        @test (mesh_u.nx, mesh_u.ny, mesh_u.nz) == (mesh_s.nx, mesh_s.ny, mesh_s.nz)
        @test maximum(maximum(abs.(mesh_u.grid[i] - mesh_s.grid[i])) for i in eachindex(mesh_s.grid)) == 0.0
        @test maximum(abs.(mesh_u.SDF .- mesh_s.SDF)) == 0.0

        # Mesh both through the shared entry point (default :linear cut points).
        generate_tetrahedral_mesh(mesh_s, joinpath(outdir, "rt_structured"))
        generate_tetrahedral_mesh(mesh_u, joinpath(outdir, "rt_unstructured"))

        # Identical topology, floating-point-close geometry (eval_sdf differs only by
        # HEX8-vs-trilinear rounding at cut points; measured node deviation ~1.4e-14).
        @test length(mesh_s.X) == BEAM_NODES && length(mesh_s.IEN) == BEAM_TETS
        @test length(mesh_u.X) == length(mesh_s.X)
        @test length(mesh_u.IEN) == length(mesh_s.IEN)
        @test all(mesh_u.IEN[i] == mesh_s.IEN[i] for i in eachindex(mesh_s.IEN))
        @test maximum(maximum(abs.(mesh_u.X[i] - mesh_s.X[i])) for i in eachindex(mesh_s.X)) < 1e-9

        # The unstructured result is independently valid.
        @test check_watertight(mesh_u).open_edges == 0
        @test count_inverted_exact(mesh_u) == 0
        @test count_components(mesh_u) == 1
        @test min_signed_volume(mesh_u) > 0
    end

    # --------------------------------------------------------------------------
    @testset "Analytic unstructured sphere" begin
        outdir = ensure_output_dir()
        center = SVector(2.0, 2.0, 2.0)
        radius = 1.2
        nodes, hexes, phi = sphere_in_distorted_block(; center = center, radius = radius)

        # Genuinely unstructured: the auto-detector must NOT mistake it for a grid.
        src = build_sdf_source(nodes, hexes, phi)
        @test src isa UnstructuredSDF
        @test bbox(src) == (SVector(0.0, 0.0, 0.0), SVector(4.0, 4.0, 4.0))

        # Source-integrated solid volume matches the analytic ball (FE quadrature).
        vol_true = 4 / 3 * pi * radius^3
        @test abs(calculate_volume_from_sdf(src) - vol_true) / vol_true <= SPHERE_U_VOL_RELERR

        # Mesh through the FE path; :bisection tracks the interpolated zero set.
        mesh = BlockMesh(src; dx = 0.15, padding = 2)
        generate_tetrahedral_mesh(mesh, joinpath(outdir, "sphere_unstructured");
                                  options = MeshGenerationOptions(cut_points = :bisection))

        # Correctness invariants (the only hard checks without an oracle).
        @test check_watertight(mesh).open_edges == 0
        @test count_inverted_exact(mesh) == 0
        @test count_components(mesh) == 1
        @test all_finite_coords(mesh)
        @test min_signed_volume(mesh) > 0

        # Quality floor (no slivers worse than the structured runs) + frozen counts.
        dih = dihedral_stats(mesh)
        @test dih.inverted == 0
        @test dih.lt5 == 0
        @test dih.min > 10.0
        @test length(mesh.X) == SPHERE_U_NODES
        @test length(mesh.IEN) == SPHERE_U_TETS

        # Surface fidelity to the TRUE sphere (the physically meaningful metric). Note
        # |eval_sdf| on the boundary is NOT ~0 here, unlike the structured paths: this
        # graded+jittered hex block is not a perfect tiling (warped, non-coplanar HEX8
        # faces make edge-adjacent elements overlap by ~1e-5), so a boundary cut point
        # can sit strictly inside two elements with different field values and eval_sdf
        # must pick one. That is a property of the distorted input, not a fidelity
        # failure -- the vertex still sits on the sphere -- so we measure GEOMETRIC
        # distance. (eval_sdf picks the most-central claimer; see UnstructuredSDF.)
        bverts = boundary_vertex_indices(mesh)
        geo = maximum(abs(norm(mesh.X[v] - center) - radius) for v in bverts)
        @test geo <= SPHERE_U_GEO_FIDELITY

        # Tet-mesh volume tracks the analytic ball.
        tetvol = sum(
            dot(mesh.X[t[2]] - mesh.X[t[1]],
                cross(mesh.X[t[3]] - mesh.X[t[1]], mesh.X[t[4]] - mesh.X[t[1]])) / 6
            for t in mesh.IEN)
        @test abs(tetvol - vol_true) / vol_true <= SPHERE_U_VOL_RELERR
    end
end
