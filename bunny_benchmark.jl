# ==============================================================================
# Stanford bunny benchmark
# ==============================================================================
#
# Meshes the Stanford bunny from its signed distance field and renders the
# dissertation figures as VECTOR PDFs (plus a PNG preview of each). They are
# meant as one row of three, in this order, sharing one camera:
#
#     output/bunny/bunny_sdf.pdf      -- 1/3 the input surface, solid grey
#     output/bunny/bunny_tet_cut.pdf  -- 2/3 cut open, ghost of the removed part
#     output/bunny/bunny_tet.pdf      -- 3/3 the whole mesh, surface with edges
#
# Run:
#     julia -t auto --project=. bunny_benchmark.jl [dx] [--no-render]
#
# `dx` is the spacing of the A15 generation lattice, i.e. THE ELEMENT SIZE, and
# it is the only knob the figures really need. The bunny is ~86 units long:
#
#     dx    cells across   nodes    tets     boundary triangles
#     8.0        11         3141    13291           2928
#     6.0        14         6853    30810           5358      <- default
#     4.0        22        20667    99933          12000
#
# All of these come out watertight, single-component, no inverted elements, with
# the smallest dihedral angle at ~14 deg. Coarser = the discretisation stays
# readable when the figure is printed small; dx = 8 is chosen for three figures
# side by side (~5 cm each), dx = 6 suits a half-page figure.
#
# ------------------------------------------------------------------------------
# Input data
# ------------------------------------------------------------------------------
# data/bunny/Z_bunny_Fine{Grid,SDF}_B-1.0.jld2, produced once from the Stanford
# bunny STL by `bunny_sdf_from_stl.jl` (see that file's header). The field is a
# NARROW BAND: stl2sdf truncates the distance at one cell size, so ~94% of the
# grid values are saturated at +-1.0042 and only the shell around the surface
# carries a true distance.
#
# That is why the SDF grid (1.0) and the meshing lattice (dx) are decoupled here:
# the source is sampled through `eval_sdf`, so the coarse lattice reads the fine
# field instead of being tied to it. It is also why `cut_points = :bisection` is
# the default -- with a saturated field the quartet-style linear estimate
# `alpha = phi_i/(phi_i - phi_j)` collapses to the edge midpoint and the surface
# would come out blocky; bisection locates the true zero of the interpolated
# field along the edge (Labelle & Shewchuk 2007, section 3.1).
# ==============================================================================

using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.Utils
using JLD2
using StaticArrays
using Printf

include(joinpath(@__DIR__, "bunny_figures.jl"))     # vector-PDF renderer

# ------------------------------------------------------------------------------
# Parameters -- mesh
# ------------------------------------------------------------------------------
const DX = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : 6.0   # element size
const CUT_POINTS = :bisection      # :bisection (see header) or :linear
const RENDER = !("--no-render" in ARGS)

const BUNNY_DIR = joinpath(@__DIR__, "data", "bunny")
const GRID_FILE = joinpath(BUNNY_DIR, "Z_bunny_FineGrid_B-1.0.jld2")
const SDF_FILE = joinpath(BUNNY_DIR, "Z_bunny_FineSDF_B-1.0.jld2")
const OUT_DIR = joinpath(@__DIR__, "output", "bunny")

# ------------------------------------------------------------------------------
# Parameters -- figures
# ------------------------------------------------------------------------------
# Source of the input-surface figure; only this one file lives outside the repo
# (the STL is git-ignored here). Missing STL = that figure is skipped.
const STL_FILE = get(ENV, "BUNNY_STL", "/Users/ondra/github/stl2sdf.jl/data/Bunny.stl")

const AZIMUTH = 45.0               # view direction, degrees (0 = looking along -Y)
const ELEVATION = 12.0             # shared by all three figures -- see below
const PAGE = 145                   # page size in points -- set it to the SIZE THE
# FIGURE IS PRINTED AT, so that LINEWIDTH (also
# in points) is not rescaled by \includegraphics.
# 145 pt ~= one of three figures across a 16 cm
# text width; use ~240 for a half-width figure.
const LINEWIDTH = 0.3              # element edges, points
const EDGE_ALPHA = 0.8             # element edges, opacity

const SOLID = RGBf(0.72, 0.72, 0.75)  # solid colour of the input surface
const FACE = RGBf(0.87, 0.87, 0.89)   # element faces -- light, so the edges read
const INNER = RGBf(0.70, 0.77, 0.92)  # faces the cut exposed -- MUST differ from
# FACE, that contrast is what makes the
# interior structure legible
const EDGE = RGBf(0.13, 0.16, 0.60)   # element edges
const GHOST_ALPHA = 0.13           # translucent shell in the cut figure

# The cut: `:crinkle` keeps whole tetrahedra (the exposed surface is made of real
# elements, stepped); `:planar` cuts the elements and gives a flat face whose
# polygons are the tetrahedra's cross-sections (calmer, shows the A15 lattice).
const CUT_MODE = :crinkle
const CUT_AZIMUTH = AZIMUTH + 35   # normal of the cutting plane, degrees
const CUT_OFFSET = 0.05            # plane shift from the centre, fraction of size

if !isfile(GRID_FILE) || !isfile(SDF_FILE)
    error("""
    Missing bunny SDF data:
        $GRID_FILE
        $SDF_FILE
    Generate it once with the sibling stl2sdf package:
        julia -t auto --project=/Users/ondra/github/stl2sdf.jl bunny_sdf_from_stl.jl
    """)
end
mkpath(OUT_DIR)

# ------------------------------------------------------------------------------
# 1. Load the field and wrap it as an SDF source
# ------------------------------------------------------------------------------
@info "Loading bunny SDF..."
@load GRID_FILE fine_grid
@load SDF_FILE fine_sdf

# The stored convention is POSITIVE = inside (as in data/beam, data/gripper);
# every SDFSource must hold phi < 0 = inside, so negate here. This is the same
# single negation `BlockMesh(fine_sdf, fine_grid)` performs internally.
grid = Array{SVector{3,Float64},3}(undef, size(fine_grid))
for i in eachindex(fine_grid)
    grid[i] = SVector{3,Float64}(fine_grid[i]...)
end
phi = -Float64.(fine_sdf)
source = StructuredSDF(grid, phi)

(bmin, bmax) = bbox(source)
@info @sprintf(
    "SDF grid %d x %d x %d, spacing %.4f, extent %.1f x %.1f x %.1f",
    size(grid)...,
    grid[2, 1, 1][1] - grid[1, 1, 1][1],
    (bmax .- bmin)...
)

# ------------------------------------------------------------------------------
# 2. Mesh it on a coarse A15 lattice
# ------------------------------------------------------------------------------
@info @sprintf(
    "Meshing with dx = %.3f (~%.0f cells across the bunny), cut_points = %s",
    DX,
    maximum(bmax .- bmin) / DX,
    CUT_POINTS
)

mesh = BlockMesh(source; dx = DX, padding = 2)
prefix = joinpath(OUT_DIR, "bunny")

t_mesh = @elapsed generate_tetrahedral_mesh(
    mesh,
    prefix;
    options = MeshGenerationOptions(scheme = "A15", cut_points = CUT_POINTS),
)

vtu_file = "$(prefix)_TriMesh-A15.vtu"

# ------------------------------------------------------------------------------
# 3. Benchmark report
# ------------------------------------------------------------------------------
include(joinpath(@__DIR__, "test", "helpers", "helpers.jl"))   # check_* / dihedral_stats

wt = check_watertight(mesh)
dh = dihedral_stats(mesh)

println()
println("=" ^ 62)
println("Stanford bunny benchmark")
println("=" ^ 62)
@printf(
    "  lattice spacing dx     %.3f  (%d x %d x %d lattice)\n",
    DX,
    mesh.nx,
    mesh.ny,
    mesh.nz
)
@printf("  cut points             %s\n", CUT_POINTS)
@printf("  nodes / tetrahedra     %d / %d\n", length(mesh.X), length(mesh.IEN))
@printf("  meshing time           %.2f s\n", t_mesh)
@printf("  open boundary edges    %d   (0 = watertight)\n", wt.open_edges)
@printf(
    "  non-manifold edges     %d   (benign thin-feature pinches)\n",
    wt.nonmanifold_edges
)
@printf("  boundary faces         %d\n", wt.boundary_faces)
@printf("  inverted elements      %d\n", count_inverted_exact(mesh))
@printf("  connected components   %d\n", count_components(mesh))
@printf("  dihedral angle range   %.2f deg .. %.2f deg\n", dh.min, dh.max)
@printf("  angles < 10 deg        %d\n", dh.lt10)
@printf("  angles > 140 deg       %d\n", dh.gt140)
println("=" ^ 62)
println("  mesh: $vtu_file")

# ------------------------------------------------------------------------------
# 4. Figures (vector PDF, + a PNG preview of each)
# ------------------------------------------------------------------------------
RENDER || exit(0)

include(joinpath(@__DIR__, "bunny_figures.jl"))

println()
@info "Rendering figures..."

draw(file, layers; lw = LINEWIDTH) = render_pdf(
    joinpath(OUT_DIR, file),
    layers;
    azimuth = AZIMUTH,
    elevation = ELEVATION,
    page = PAGE,
    linewidth = lw,
)

# -- 1/3: the input surface -- the Stanford bunny STL the SDF was sampled from --
if isfile(STL_FILE)
    draw(
        "bunny_sdf.pdf",
        [Layer(faces = read_ascii_stl(STL_FILE), color = SOLID)];
        lw = 0.25,
    )
else
    @warn "STL not found, skipping the input-surface figure: $STL_FILE"
end

# -- 2/3: cut open, removed part left as a translucent ghost -------------------
# The cutting plane is turned relative to the camera rather than the object: all
# three figures of the row must share one view, or the row reads as a jumble.
cut_normal = SVector(sind(CUT_AZIMUTH), -cosd(CUT_AZIMUTH), 0.0)
cut_origin = (bmin + bmax) / 2 + CUT_OFFSET * maximum(bmax - bmin) * cut_normal
removed_side(p) = dot(p - cut_origin, cut_normal) > 0

outer_keys, surface = boundary_faces(mesh.X, mesh.IEN)
kept, _ = clip_tets(mesh.X, mesh.IEN, cut_origin, cut_normal)

if CUT_MODE === :crinkle
    # Whole elements: the cut keeps every tet whose centroid stays on the near
    # side, so the exposed surface is made of real tetrahedra.
    skin, exposed = split_boundary(mesh.X, kept, Set(outer_keys))
else
    # True planar section: elements are cut, and each polygon of the flat face is
    # one tetrahedron's cross-section.
    skin = filter(!isnothing, [clip_face(f, cut_origin, cut_normal) for f in surface])
    exposed = section_faces(mesh.X, mesh.IEN, cut_origin, cut_normal)
end

# The ghost is the OUTER skin of the removed part, taken from the full mesh's
# boundary, so it carries no cut face of its own to fight with the opaque one.
ghost = filter(f -> removed_side(sum(f) / length(f)), surface)

draw(
    "bunny_tet_cut.pdf",
    [
        Layer(faces = skin, color = FACE, stroke = EDGE, stroke_alpha = EDGE_ALPHA),
        Layer(faces = exposed, color = INNER, stroke = EDGE, stroke_alpha = EDGE_ALPHA),
        Layer(faces = ghost, color = FACE, alpha = GHOST_ALPHA, stroke = EDGE),
    ],
)

# -- 3/3: the whole mesh, surface with edges -----------------------------------
draw(
    "bunny_tet.pdf",
    [Layer(faces = surface, color = FACE, stroke = EDGE, stroke_alpha = EDGE_ALPHA)],
)

println("\nFigures in $OUT_DIR:")
for f in sort(filter(f -> endswith(f, ".pdf"), readdir(OUT_DIR)))
    @printf("  %-22s %6.2f MB\n", f, filesize(joinpath(OUT_DIR, f)) / 2^20)
end
