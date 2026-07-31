# ==============================================================================
# Stanford bunny: STL -> signed distance field  (one-off data preparation)
# ==============================================================================
#
# Produces the two JLD2 files that `bunny_benchmark.jl` consumes:
#
#     data/bunny/Z_bunny_FineGrid_B-<B>.jld2   ->  fine_grid
#     data/bunny/Z_bunny_FineSDF_B-<B>.jld2    ->  fine_sdf
#
# stored in exactly the same layout and sign convention as the beam / gripper
# data already in `data/` (POSITIVE = inside; `BlockMesh` negates on load).
#
# This step needs the sibling package stl2sdf.jl, which is NOT a dependency of
# Implicit2TetMesh -- run it in that project's environment:
#
#     julia -t auto --project=/Users/ondra/github/stl2sdf.jl bunny_sdf_from_stl.jl
#
# Optional positional arguments:  <cell_size>  <stl_file>
#
# It only has to be run once; the resulting field is resolution-independent
# input for the benchmark, which picks its own (coarser) meshing lattice.
# ==============================================================================

using stl2sdf
using stl2sdf.DataImport
using stl2sdf.MeshGrid
using stl2sdf.SignedDistances
using StaticArrays
using JLD2

# ------------------------------------------------------------------------------
# Parameters
# ------------------------------------------------------------------------------
const CELL_SIZE = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 1.0
const STL_FILE = length(ARGS) >= 2 ? ARGS[2] :
                 "/Users/ondra/github/stl2sdf.jl/data/Bunny.stl"
const OUT_DIR = joinpath(@__DIR__, "data", "bunny")

isfile(STL_FILE) || error("STL not found: $STL_FILE")
mkpath(OUT_DIR)

println("STL:        $STL_FILE")
println("cell size:  $CELL_SIZE")
println("threads:    $(Threads.nthreads())")

# ------------------------------------------------------------------------------
# 1. Import the triangulated surface
# ------------------------------------------------------------------------------
(X, IEN) = import_stl(STL_FILE)
tri_mesh = TriangularMesh(X, IEN)
println("surface:    $(tri_mesh.nnp) vertices, $(tri_mesh.nel) triangles")

# ------------------------------------------------------------------------------
# 2. Sampling grid (uniform, 3 margin cells around the surface AABB)
# ------------------------------------------------------------------------------
sdf_grid = noninteractive_sdf_grid_setup(tri_mesh, CELL_SIZE)
points = generateGridPoints(sdf_grid)

# ------------------------------------------------------------------------------
# 3. Unsigned distance + sign (ray casting with winding-number fallback)
# ------------------------------------------------------------------------------
t_dist = @elapsed ((dists, _xp) = evalDistancesOnTriMesh(tri_mesh, sdf_grid, points))
println("distances:  $(round(t_dist, digits = 1)) s")

t_sign = @elapsed ((signs, confidences) = raycast_sign_detection(tri_mesh, sdf_grid, points))
println("signs:      $(round(t_sign, digits = 1)) s")
low = count(c -> c < 0.6, confidences)
low == 0 || @warn "$low grid points with low sign confidence (<0.6)"

# POSITIVE = inside, matching data/beam and data/gripper.
sdf_values = dists .* signs

# ------------------------------------------------------------------------------
# 4. Reshape to the structured layout used by `BlockMesh(fine_sdf, fine_grid)`
# ------------------------------------------------------------------------------
# `evalDistancesOnTriMesh` leaves an unreached point at the sentinel -1e10; clamp
# it to the largest real magnitude, as stl2sdf's own export does.
sdf_f32 = Float32.(sdf_values)
finite_max = maximum(abs, filter(x -> abs(x) < 1.0f9, sdf_f32))
@. sdf_f32 = ifelse(abs(sdf_f32) ≈ 1.0f10, sign(sdf_f32) * finite_max, sdf_f32)

dims = Tuple(sdf_grid.N .+ 1)                     # grid points per axis (i fastest)
fine_sdf = reshape(sdf_f32, dims)

xs = range(Float32(sdf_grid.AABB_min[1]), Float32(sdf_grid.AABB_max[1]), length = dims[1])
ys = range(Float32(sdf_grid.AABB_min[2]), Float32(sdf_grid.AABB_max[2]), length = dims[2])
zs = range(Float32(sdf_grid.AABB_min[3]), Float32(sdf_grid.AABB_max[3]), length = dims[3])
# SVector (not Vector) keeps the file ~10x smaller; BlockMesh splats either one.
fine_grid = [SVector{3,Float32}(x, y, z) for x in xs, y in ys, z in zs]

# ------------------------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------------------------
B = round(CELL_SIZE, digits = 4)
grid_file = joinpath(OUT_DIR, "Z_bunny_FineGrid_B-$(B).jld2")
sdf_file = joinpath(OUT_DIR, "Z_bunny_FineSDF_B-$(B).jld2")

@save grid_file fine_grid
@save sdf_file fine_sdf

println("\ngrid:       $(dims) points, spacing $(round(sdf_grid.cell_size, digits = 4))")
println("sdf range:  $(extrema(fine_sdf))  (positive = inside)")
println("written:    $grid_file")
println("            $sdf_file")
