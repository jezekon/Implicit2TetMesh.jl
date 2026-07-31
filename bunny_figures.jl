# ==============================================================================
# Vector-PDF figures for the Stanford bunny benchmark
# ==============================================================================
#
# Included by `bunny_benchmark.jl`; renders a `BlockMesh` (and, optionally, the
# input STL surface) straight to PDF with CairoMakie.
#
# WHY NOT PARAVIEW: ParaView's "Export Scene -> PDF/SVG" runs through GL2PS, and
# on the OpenGL2 backend it embeds a raster image of the viewport no matter what
# `Rasterize 3D geometry` is set to -- the exported PDF is a bitmap in a PDF
# wrapper (verified on ParaView 6.1.1 / macOS). So the figures are drawn here
# instead: triangles are projected, flat-shaded, sorted back-to-front and emitted
# as filled paths, which is genuine vector output and gives exact control over
# element-edge width -- the whole point of the coarse-mesh figure.
#
# Hidden-surface removal is the painter's algorithm (per-triangle fill + stroke,
# far to near), which is exact for a closed surface whose triangles do not
# intersect -- true for a boundary triangulation.
# ==============================================================================

using CairoMakie
using CairoMakie.GeometryBasics
using LinearAlgebra
using StaticArrays

const Tri = NTuple{3,SVector{3,Float64}}

# ------------------------------------------------------------------------------
# Geometry extraction
# ------------------------------------------------------------------------------

"""
    boundary_triangles(X, tets) -> Vector{Tri}

Outward-oriented boundary triangles of a tetrahedral element list: every face
shared by exactly one tet, wound so that its normal points away from the tet's
fourth vertex.
"""
function boundary_triangles(X::Vector{SVector{3,Float64}}, tets)
    count = Dict{NTuple{3,Int},Int}()
    owner = Dict{NTuple{3,Int},NTuple{4,Int}}()
    for tet in tets
        a, b, c, d = tet[1], tet[2], tet[3], tet[4]
        for (i, j, k, l) in ((a, b, c, d), (a, b, d, c), (a, c, d, b), (b, c, d, a))
            key = extrema3(i, j, k)
            count[key] = get(count, key, 0) + 1
            owner[key] = (i, j, k, l)
        end
    end

    tris = Tri[]
    for (key, c) in count
        c == 1 || continue
        (i, j, k, l) = owner[key]
        p, q, r = X[i], X[j], X[k]
        n = cross(q - p, r - p)
        # Flip the winding if the normal points at the opposite vertex (inward).
        dot(n, X[l] - p) > 0 ? push!(tris, (p, r, q)) : push!(tris, (p, q, r))
    end
    return tris
end

extrema3(a::Int, b::Int, c::Int) =
    (min(a, b, c), a + b + c - min(a, b, c) - max(a, b, c), max(a, b, c))

"""
    clip_tets(X, tets, origin, normal) -> kept, removed

Split an element list by a plane WITHOUT cutting elements ("crinkle clip"): a tet
goes to `kept` when its centroid is on the negative side of the plane.
"""
function clip_tets(X::Vector{SVector{3,Float64}}, tets, origin::SVector{3,Float64},
                   normal::SVector{3,Float64})
    kept = eltype(tets)[]
    removed = eltype(tets)[]
    for tet in tets
        ctr = (X[tet[1]] + X[tet[2]] + X[tet[3]] + X[tet[4]]) / 4
        push!(dot(ctr - origin, normal) <= 0 ? kept : removed, tet)
    end
    return kept, removed
end

"""
    read_ascii_stl(path) -> Vector{Tri}

Minimal reader for the ASCII STL the bunny ships as -- vertex lines only, the
facet normals are recomputed from the winding when shading.
"""
function read_ascii_stl(path::String)
    tris = Tri[]
    buf = SVector{3,Float64}[]
    for line in eachline(path)
        s = lstrip(line)
        startswith(s, "vertex") || continue
        f = split(s)
        push!(buf, SVector{3,Float64}(parse(Float64, f[2]), parse(Float64, f[3]),
                                      parse(Float64, f[4])))
        if length(buf) == 3
            push!(tris, (buf[1], buf[2], buf[3]))
            empty!(buf)
        end
    end
    return tris
end

# ------------------------------------------------------------------------------
# Rendering
# ------------------------------------------------------------------------------

"""
    Layer(tris; color, alpha, stroke, cull)

One group of triangles with its own style. `stroke = nothing` draws no element
edges; `cull = false` keeps back faces (wanted for a translucent shell).
"""
Base.@kwdef struct Layer
    tris::Vector{Tri}
    color::RGBf = RGBf(0.13, 0.16, 0.60)
    alpha::Float64 = 1.0
    stroke::Union{Nothing,RGBf} = nothing
    cull::Bool = true
end

"""
    render_pdf(file, layers; kwargs...)

Project, shade, depth-sort and draw `layers`, then write `file` (.pdf keeps the
output vector; .png is written too when `preview = true`).

Keyword arguments:
- `azimuth`, `elevation`: view direction in degrees; azimuth 0 looks along -Y,
  positive azimuth turns around +Z (up).
- `page`: page size in PostScript points. Make it the size the figure is printed
  at, so `linewidth` (also in points) means what it says.
- `linewidth`: element-edge width in points.
- `zoom`: >1 crops closer.
"""
function render_pdf(file::String, layers::Vector{Layer};
                    azimuth::Real = 235.0, elevation::Real = 14.0,
                    page::Real = 240, linewidth::Real = 0.3, zoom::Real = 1.0,
                    ambient::Real = 0.42, diffuse::Real = 0.58, preview::Bool = true)

    # Camera frame: `dir` points from the scene towards the viewer.
    a, e = deg2rad(azimuth), deg2rad(elevation)
    dir = SVector(sin(a) * cos(e), -cos(a) * cos(e), sin(e))
    right = normalize(cross(SVector(0.0, 0.0, 1.0), dir))
    up = cross(dir, right)
    light = normalize(dir + 0.45 * right + 0.35 * up)

    polys = Polygon{2,Float32}[]
    fills = RGBAf[]
    strokes = RGBAf[]
    depths = Float64[]

    for layer in layers
        base = layer.color
        for t in layer.tris
            n = cross(t[2] - t[1], t[3] - t[1])
            nn = norm(n)
            nn > 1e-14 || continue
            n /= nn
            facing = dot(n, dir)
            (layer.cull && facing <= 0) && continue
            # Two-sided shading, so back faces of a translucent shell are lit too.
            s = clamp(ambient + diffuse * abs(dot(n, light)), 0.0, 1.0)

            pts = ntuple(i -> Point2f(dot(t[i], right), dot(t[i], up)), 3)
            fill = RGBAf(base.r * s, base.g * s, base.b * s, layer.alpha)
            push!(polys, Polygon(collect(pts)))
            push!(fills, fill)
            # A stroke is always drawn: in its own colour it marks element edges,
            # in the fill colour it closes the hairline seams Cairo leaves between
            # adjacent filled paths (which otherwise texture a smooth surface).
            push!(strokes, layer.stroke === nothing ? fill :
                           RGBAf(layer.stroke.r, layer.stroke.g, layer.stroke.b, layer.alpha))
            push!(depths, (dot(t[1], dir) + dot(t[2], dir) + dot(t[3], dir)) / 3)
        end
    end
    isempty(polys) && error("render_pdf: nothing visible to draw")

    order = sortperm(depths)                 # far first -- painter's algorithm

    fig = Figure(size = (page, page), backgroundcolor = :white, figure_padding = 0)
    ax = Axis(fig[1, 1], aspect = DataAspect(), backgroundcolor = :white)
    hidedecorations!(ax)
    hidespines!(ax)
    ax.xautolimitmargin = (0.02, 0.02)
    ax.yautolimitmargin = (0.02, 0.02)

    poly!(ax, polys[order]; color = fills[order], strokecolor = strokes[order],
          strokewidth = linewidth)

    if zoom != 1
        autolimits!(ax)
        lims = ax.finallimits[]
        cx, cy = lims.origin[1] + lims.widths[1] / 2, lims.origin[2] + lims.widths[2] / 2
        w, h = lims.widths[1] / zoom, lims.widths[2] / zoom
        limits!(ax, cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2)
    end

    save(file, fig)
    preview && save(replace(file, r"\.pdf$" => ".png"), fig, px_per_unit = 4)
    return length(polys)
end
