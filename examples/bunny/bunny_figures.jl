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
# instead: faces are projected, flat-shaded, sorted back-to-front and emitted as
# filled paths, which is genuine vector output and gives exact control over
# element-edge width -- the whole point of the coarse-mesh figure.
#
# Hidden-surface removal is the painter's algorithm (per-face fill + stroke, far
# to near), which is exact for surfaces whose faces do not intersect -- true for
# a boundary triangulation and for a planar section of it.
# ==============================================================================

using CairoMakie
using CairoMakie.GeometryBasics
using LinearAlgebra
using StaticArrays

# A face is a planar convex polygon: 3 points for a boundary triangle, 3 or 4 for
# a tetrahedron's planar cross-section. Sections are kept as polygons rather than
# fan-triangulated so that stroking them does not draw a spurious diagonal.
const Face = Vector{SVector{3,Float64}}

# ------------------------------------------------------------------------------
# Boundary extraction
# ------------------------------------------------------------------------------

facekey(a::Int, b::Int, c::Int) =
    (min(a, b, c), a + b + c - min(a, b, c) - max(a, b, c), max(a, b, c))

"""
    boundary_faces(X, tets) -> (keys, faces)

Outward-oriented boundary triangles of a tetrahedral element list -- every face
used by exactly one tet, wound so its normal points away from that tet's fourth
vertex -- together with the sorted vertex-index key of each, so that two element
sets can be compared face by face.
"""
function boundary_faces(X::Vector{SVector{3,Float64}}, tets)
    count = Dict{NTuple{3,Int},Int}()
    owner = Dict{NTuple{3,Int},NTuple{4,Int}}()
    for tet in tets
        a, b, c, d = tet[1], tet[2], tet[3], tet[4]
        for (i, j, k, l) in ((a, b, c, d), (a, b, d, c), (a, c, d, b), (b, c, d, a))
            key = facekey(i, j, k)
            count[key] = get(count, key, 0) + 1
            owner[key] = (i, j, k, l)
        end
    end

    keys = NTuple{3,Int}[]
    faces = Face[]
    for (key, c) in count
        c == 1 || continue
        (i, j, k, l) = owner[key]
        p, q, r = X[i], X[j], X[k]
        inward = dot(cross(q - p, r - p), X[l] - p) > 0
        push!(keys, key)
        push!(faces, inward ? Face([p, r, q]) : Face([p, q, r]))
    end
    return keys, faces
end

boundary_triangles(X, tets) = boundary_faces(X, tets)[2]

"""
    split_boundary(X, tets, outer_keys) -> (outer, exposed)

Boundary faces of `tets` split into those that are also faces of the full mesh
(`outer_keys`, the object's own skin) and those the cut has newly exposed. Giving
the two groups different colours is what makes the interior structure read.
"""
function split_boundary(X, tets, outer_keys::Set{NTuple{3,Int}})
    keys, faces = boundary_faces(X, tets)
    outer = Face[]
    exposed = Face[]
    for (k, f) in zip(keys, faces)
        push!(k in outer_keys ? outer : exposed, f)
    end
    return outer, exposed
end

# ------------------------------------------------------------------------------
# Cutting
# ------------------------------------------------------------------------------

"""
    clip_tets(X, tets, origin, normal) -> kept, removed

Split an element list by a plane WITHOUT cutting elements ("crinkle clip"): a tet
is kept when its centroid is on the negative side of the plane.
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
    clip_tets(X, tets, keep) -> kept, removed

Crinkle clip against an arbitrary region: `keep(centroid)` decides. Use this for
anything a single plane cannot express, e.g. a wedge cut out of a corner.
"""
function clip_tets(X::Vector{SVector{3,Float64}}, tets, keep::Function)
    kept = eltype(tets)[]
    removed = eltype(tets)[]
    for tet in tets
        ctr = (X[tet[1]] + X[tet[2]] + X[tet[3]] + X[tet[4]]) / 4
        push!(keep(ctr) ? kept : removed, tet)
    end
    return kept, removed
end

"""
    clip_face(f, origin, normal) -> Union{Face,Nothing}

Sutherland-Hodgman clip of a convex face against the half-space
`dot(p - origin, normal) <= 0`. `nothing` when nothing survives.
"""
function clip_face(f::Face, origin::SVector{3,Float64}, normal::SVector{3,Float64})
    out = SVector{3,Float64}[]
    n = length(f)
    for i = 1:n
        p, q = f[i], f[mod1(i + 1, n)]
        dp, dq = dot(p - origin, normal), dot(q - origin, normal)
        dp <= 0 && push!(out, p)
        if (dp < 0) != (dq < 0)
            push!(out, p + (dp / (dp - dq)) * (q - p))
        end
    end
    return length(out) >= 3 ? out : nothing
end

"""
    section_faces(X, tets, origin, normal) -> Vector{Face}

Cross-section of every tetrahedron that the plane passes through, as convex
polygons wound so their normal is `+normal` (i.e. facing out of the kept, negative
side). This is the flat cut face of a true planar clip; each polygon is one
element, so the section shows the actual internal element pattern.
"""
function section_faces(X::Vector{SVector{3,Float64}}, tets, origin::SVector{3,Float64},
                       normal::SVector{3,Float64})
    edges = ((1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4))
    faces = Face[]
    pts = SVector{3,Float64}[]
    for tet in tets
        d = ntuple(i -> dot(X[tet[i]] - origin, normal), 4)
        (minimum(d) < 0 && maximum(d) > 0) || continue
        empty!(pts)
        for (a, b) in edges
            (d[a] < 0) == (d[b] < 0) && continue
            p, q = X[tet[a]], X[tet[b]]
            push!(pts, p + (d[a] / (d[a] - d[b])) * (q - p))
        end
        length(pts) >= 3 || continue

        # Order the points around their centroid in the cut plane, then orient.
        c = sum(pts) / length(pts)
        u = normalize(pts[1] - c)
        v = cross(normal, u)
        order = sortperm([atan(dot(p - c, v), dot(p - c, u)) for p in pts])
        poly = pts[order]
        dot(cross(poly[2] - poly[1], poly[3] - poly[1]), normal) < 0 && reverse!(poly)
        push!(faces, poly)
    end
    return faces
end

"""
    read_ascii_stl(path) -> Vector{Face}

Minimal reader for the ASCII STL the bunny ships as -- vertex lines only; facet
normals are recomputed from the winding when shading.
"""
function read_ascii_stl(path::String)
    faces = Face[]
    buf = SVector{3,Float64}[]
    for line in eachline(path)
        s = lstrip(line)
        startswith(s, "vertex") || continue
        f = split(s)
        push!(buf, SVector{3,Float64}(parse(Float64, f[2]), parse(Float64, f[3]),
                                      parse(Float64, f[4])))
        if length(buf) == 3
            push!(faces, copy(buf))
            empty!(buf)
        end
    end
    return faces
end

# ------------------------------------------------------------------------------
# Rendering
# ------------------------------------------------------------------------------

"""
    Layer(faces; color, alpha, stroke, stroke_alpha, cull)

One group of faces with its own style. `stroke = nothing` draws no element edges;
`cull = false` keeps back faces (wanted for a translucent shell).
"""
Base.@kwdef struct Layer
    faces::Vector{Face}
    color::RGBf = RGBf(0.87, 0.87, 0.89)
    alpha::Float64 = 1.0
    stroke::Union{Nothing,RGBf} = nothing
    stroke_alpha::Float64 = 1.0
    cull::Bool = true
end

"""
    render_pdf(file, layers; kwargs...)

Project, shade, depth-sort and draw `layers`, then write `file` (.pdf stays
vector; a .png preview is written alongside when `preview = true`).

Keyword arguments:
- `azimuth`, `elevation`: view direction in degrees; azimuth 0 looks along -Y,
  positive azimuth turns around +Z (up).
- `page`: page size in PostScript points. Make it the size the figure is printed
  at, so `linewidth` (also in points) means what it says.
- `linewidth`: element-edge width in points.
- `zoom`: >1 crops closer.
"""
function render_pdf(file::String, layers::Vector{Layer};
                    azimuth::Real = 45.0, elevation::Real = 12.0,
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
        for f in layer.faces
            n = cross(f[2] - f[1], f[3] - f[1])
            nn = norm(n)
            nn > 1e-14 || continue
            n /= nn
            (layer.cull && dot(n, dir) <= 0) && continue
            # Two-sided shading, so back faces of a translucent shell are lit too.
            s = clamp(ambient + diffuse * abs(dot(n, light)), 0.0, 1.0)

            fill = RGBAf(base.r * s, base.g * s, base.b * s, layer.alpha)
            push!(polys, Polygon([Point2f(dot(p, right), dot(p, up)) for p in f]))
            push!(fills, fill)
            # A stroke is always drawn: in its own colour it marks element edges,
            # in the fill colour it closes the hairline seams Cairo leaves between
            # adjacent filled paths (which otherwise texture a smooth surface).
            push!(strokes, layer.stroke === nothing ? fill :
                  RGBAf(layer.stroke.r, layer.stroke.g, layer.stroke.b,
                        layer.alpha * layer.stroke_alpha))
            push!(depths, sum(dot(p, dir) for p in f) / length(f))
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
