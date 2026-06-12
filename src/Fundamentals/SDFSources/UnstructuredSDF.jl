# ----------------------------------------------------------------------------
# Unstructured HEX8 SDF source (finite-element shape-function interpolation)
# ----------------------------------------------------------------------------
# Evaluates an implicit field defined at the nodes of a conforming 8-node
# hexahedral (HEX8) mesh, exactly as a trilinear FE field: locate the element
# containing the query point, invert its isoparametric map to natural coordinates
# (xi, eta, zeta) in [-1, 1]^3, then interpolate the nodal values with the HEX8
# shape functions. On a CONFORMING mesh the restriction of this field to a shared
# face depends only on the four face nodes, so the field is C0 and the zero
# isosurface is crack-free -- the same property the structured trilinear path has.

# Natural (reference) coordinates of the 8 HEX8 corner nodes, in the standard VTK
# ordering. This is the SAME corner ordering get_cell_sdf_values uses for a lattice
# cell, so a structured grid re-expressed as HEX8 elements maps consistently.
const HEX8_NATURAL = (
    SVector(-1.0, -1.0, -1.0),   # 1
    SVector(1.0, -1.0, -1.0),    # 2
    SVector(1.0, 1.0, -1.0),     # 3
    SVector(-1.0, 1.0, -1.0),    # 4
    SVector(-1.0, -1.0, 1.0),    # 5
    SVector(1.0, -1.0, 1.0),     # 6
    SVector(1.0, 1.0, 1.0),      # 7
    SVector(-1.0, 1.0, 1.0),     # 8
)

"""
    hex8_shape(xi) -> SVector{8,Float64}

Trilinear HEX8 shape functions at natural coordinates `xi = (xi, eta, zeta)`,
`N_a = 1/8 (1 + xi*xi_a)(1 + eta*eta_a)(1 + zeta*zeta_a)`, in VTK node order. Used
to interpolate the nodal field (and node coordinates) inside an element.
"""
function hex8_shape(xi::SVector{3,Float64})
    x, y, z = xi[1], xi[2], xi[3]
    return SVector{8,Float64}(
        0.125 * (1 - x) * (1 - y) * (1 - z),
        0.125 * (1 + x) * (1 - y) * (1 - z),
        0.125 * (1 + x) * (1 + y) * (1 - z),
        0.125 * (1 - x) * (1 + y) * (1 - z),
        0.125 * (1 - x) * (1 - y) * (1 + z),
        0.125 * (1 + x) * (1 - y) * (1 + z),
        0.125 * (1 + x) * (1 + y) * (1 + z),
        0.125 * (1 - x) * (1 + y) * (1 + z),
    )
end

"""
    hex8_map_and_jacobian(X, xi) -> (point, J)

Physical point `sum_a N_a(xi) X_a` and the 3x3 Jacobian `J[m,n] = d(point_m)/d(xi_n)`
of the isoparametric map for a HEX8 element with corner coordinates `X` (VTK order)
at natural coordinates `xi`. The Jacobian drives the Newton inverse map and the
volume quadrature.
"""
function hex8_map_and_jacobian(X::NTuple{8,SVector{3,Float64}}, xi::SVector{3,Float64})
    point = zero(SVector{3,Float64})
    J = zero(SMatrix{3,3,Float64})
    for a = 1:8
        na = HEX8_NATURAL[a]
        sa = 1 + xi[1] * na[1]
        ta = 1 + xi[2] * na[2]
        ua = 1 + xi[3] * na[3]
        N = 0.125 * sa * ta * ua
        # Shape-function gradient w.r.t. the natural coordinates.
        dN = SVector(0.125 * na[1] * ta * ua, 0.125 * na[2] * sa * ua, 0.125 * na[3] * sa * ta)
        point = point + N * X[a]
        J = J + X[a] * dN'          # outer product accumulates J column by column
    end
    return point, J
end

"""
    inverse_hex_map(X, p; maxiter, xi_tol, inside_tol) -> (xi, inside)

Invert the HEX8 isoparametric map: find natural coordinates `xi` with
`sum_a N_a(xi) X_a == p` by Newton's method (`xi <- xi - J\\(map(xi) - p)`), starting
at the element center. The map is nonlinear for a general (non-parallelepiped) hex,
so a few iterations are needed. `inside` is true when the solution lies within the
reference cube `[-1,1]^3` (with tolerance `inside_tol`), i.e. `p` is inside this
element. A non-finite Newton step (degenerate/near-singular Jacobian away from the
cell) returns `inside = false`.
"""
function inverse_hex_map(
    X::NTuple{8,SVector{3,Float64}},
    p::SVector{3,Float64};
    maxiter::Int = 25,
    xi_tol::Float64 = 1e-12,
    inside_tol::Float64 = 1e-9,
)
    xi = zero(SVector{3,Float64})
    for _ = 1:maxiter
        point, J = hex8_map_and_jacobian(X, xi)
        dxi = J \ (point - p)
        xi = xi - dxi
        all(isfinite, xi) || return xi, false
        maximum(abs.(dxi)) <= xi_tol && break
    end
    inside = abs(xi[1]) <= 1 + inside_tol &&
             abs(xi[2]) <= 1 + inside_tol &&
             abs(xi[3]) <= 1 + inside_tol
    return xi, inside
end

"""
    hex_volume(X) -> Float64

Volume of a HEX8 element with corner coordinates `X`, integrated with 2x2x2
Gauss-Legendre quadrature of `|det J|` (exact for the trilinear map).
"""
function hex_volume(X::NTuple{8,SVector{3,Float64}})
    g = 1 / sqrt(3.0)
    pts = (-g, g)
    v = 0.0
    for c in pts, b in pts, a in pts            # all weights are 1
        _, J = hex8_map_and_jacobian(X, SVector(a, b, c))
        v += abs(det(J))
    end
    return v
end

# ----------------------------------------------------------------------------
# The source type + uniform spatial hash for point location
# ----------------------------------------------------------------------------

"""
    UnstructuredSDF(nodes, hexes, phi; validate = true)

Implicit field on a conforming HEX8 mesh: `nodes` are node coordinates, `hexes` the
8-index element connectivity (VTK corner order), and `phi` the nodal field already
in the phi<0=inside convention. Between nodes the field is the HEX8 shape-function
interpolation (see [`hex8_shape`](@ref)).

Point location uses a uniform spatial hash over element bounding boxes: the field
query hashes the point to a cell and tests only the few elements registered there
with the Newton inverse map. A point claimed by more than one element (a shared
face/edge, or a slight overlap between distorted hexes whose warped faces do not tile
space exactly) is evaluated in the element that contains it MOST CENTRALLY, so the
field is single-valued and robust to mildly non-tiling inputs (see [`eval_sdf`](@ref)).
The mesh MUST be conforming (no hanging nodes) and have positive Jacobians; `validate`
checks the Jacobian sign of every element at its center and the index ranges (it does
not detect hanging nodes -- conformance is the caller's contract).

Outside-domain rule: a point outside the meshed region returns a positive (exterior)
value equal to its distance to the mesh bounding box. This is exact for a box-shaped
domain -- the SIMP / design-domain case, where only the padding ring of the
generation lattice falls outside -- and keeps the zero isosurface enclosed.
"""
struct UnstructuredSDF <: SDFSource
    nodes::Vector{SVector{3,Float64}}
    hexes::Vector{NTuple{8,Int}}
    phi::Vector{Float64}
    bmin::SVector{3,Float64}             # node bounding box (also the domain box)
    bmax::SVector{3,Float64}
    nc::NTuple{3,Int}                    # spatial-hash cell counts per axis
    inv_cell::SVector{3,Float64}         # nc / extent per axis (0 on a flat axis)
    buckets::Dict{NTuple{3,Int},Vector{Int}}   # hash cell -> element indices overlapping it
end

# Map a point to its spatial-hash cell (clamped into range; a flat axis -> 0).
function hash_cell(bmin::SVector{3,Float64}, inv_cell::SVector{3,Float64},
                   nc::NTuple{3,Int}, p::SVector{3,Float64})
    ix = clamp(floor(Int, (p[1] - bmin[1]) * inv_cell[1]), 0, nc[1] - 1)
    iy = clamp(floor(Int, (p[2] - bmin[2]) * inv_cell[2]), 0, nc[2] - 1)
    iz = clamp(floor(Int, (p[3] - bmin[3]) * inv_cell[3]), 0, nc[3] - 1)
    return (ix, iy, iz)
end

# The 8 corner coordinates of element e, in VTK order.
hex_nodes(s::UnstructuredSDF, e::Int) = ntuple(a -> s.nodes[s.hexes[e][a]], 8)

# Build the uniform spatial hash: choose a cell size ~ the mean element size and
# register each element into every hash cell its bounding box overlaps.
function build_spatial_hash(nodes::Vector{SVector{3,Float64}}, hexes::Vector{NTuple{8,Int}})
    bmin = nodes[1]
    bmax = nodes[1]
    for p in nodes
        bmin = min.(bmin, p)
        bmax = max.(bmax, p)
    end

    # Characteristic element size = mean of the per-element bounding-box max extent.
    n = length(hexes)
    size_sum = 0.0
    for hx in hexes
        lo = nodes[hx[1]]
        hi = nodes[hx[1]]
        for a = 2:8
            lo = min.(lo, nodes[hx[a]])
            hi = max.(hi, nodes[hx[a]])
        end
        size_sum += maximum(hi - lo)
    end
    h = size_sum / max(n, 1)
    h > 0 || error("UnstructuredSDF: degenerate mesh (zero characteristic element size)")

    ext = bmax - bmin
    nc = ntuple(d -> clamp(round(Int, ext[d] / h), 1, 256), 3)
    inv_cell = SVector{3,Float64}(ntuple(d -> ext[d] > 0 ? nc[d] / ext[d] : 0.0, 3))

    buckets = Dict{NTuple{3,Int},Vector{Int}}()
    for (e, hx) in enumerate(hexes)
        lo = nodes[hx[1]]
        hi = nodes[hx[1]]
        for a = 2:8
            lo = min.(lo, nodes[hx[a]])
            hi = max.(hi, nodes[hx[a]])
        end
        lo_c = hash_cell(bmin, inv_cell, nc, lo)
        hi_c = hash_cell(bmin, inv_cell, nc, hi)
        for cz = lo_c[3]:hi_c[3], cy = lo_c[2]:hi_c[2], cx = lo_c[1]:hi_c[1]
            push!(get!(buckets, (cx, cy, cz), Int[]), e)
        end
    end
    return bmin, bmax, nc, inv_cell, buckets
end

"""
    UnstructuredSDF(nodes, hexes, phi; validate = true)

Convenience constructor accepting `nodes` as any vector of 3-component coordinates
and `hexes` as any vector of 8 integer indices; it normalizes them to the stored
types, validates, builds the spatial hash, and returns the source.
"""
function UnstructuredSDF(
    nodes::AbstractVector,
    hexes::AbstractVector,
    phi::AbstractVector;
    validate::Bool = true,
)
    nodes_sv = SVector{3,Float64}[SVector{3,Float64}(p[1], p[2], p[3]) for p in nodes]
    hexes_nt = NTuple{8,Int}[
        (Int(h[1]), Int(h[2]), Int(h[3]), Int(h[4]), Int(h[5]), Int(h[6]), Int(h[7]), Int(h[8])) for h in hexes
    ]
    phi_f = Float64.(phi)

    length(phi_f) == length(nodes_sv) ||
        error("UnstructuredSDF: phi length $(length(phi_f)) != node count $(length(nodes_sv))")
    isempty(hexes_nt) && error("UnstructuredSDF: empty element list")

    if validate
        nnodes = length(nodes_sv)
        for (e, hx) in enumerate(hexes_nt)
            for a = 1:8
                1 <= hx[a] <= nnodes ||
                    error("UnstructuredSDF: element $e references node $(hx[a]) outside 1:$nnodes")
            end
            X = ntuple(a -> nodes_sv[hx[a]], 8)
            _, J = hex8_map_and_jacobian(X, zero(SVector{3,Float64}))
            det(J) > 0 || error(
                "UnstructuredSDF: element $e has non-positive Jacobian (det = $(det(J))) at its " *
                "center; check the HEX8 node ordering (VTK corner order expected)",
            )
        end
    end

    bmin, bmax, nc, inv_cell, buckets = build_spatial_hash(nodes_sv, hexes_nt)
    return UnstructuredSDF(nodes_sv, hexes_nt, phi_f, bmin, bmax, nc, inv_cell, buckets)
end

"""
    bbox(s::UnstructuredSDF) -> (min_corner, max_corner)

Bounding box of the HEX8 node set (the domain box for a box-shaped mesh).
"""
bbox(s::UnstructuredSDF) = (s.bmin, s.bmax)

"""
    eval_sdf(s::UnstructuredSDF, p) -> Float64

Field value at `p`: locate the containing element via the spatial hash + Newton
inverse map and interpolate the nodal values with the HEX8 shape functions; if `p`
lies outside every element, return the (positive) distance to the domain bounding box.
See the type docstring for the outside-domain rule and its box-domain assumption.

When several elements claim `p` -- it sits on a shared face/edge, or two warped hexes
mildly overlap (non-coplanar faces of a distorted hex mesh do not tile space exactly)
-- the field is evaluated in the element that contains `p` MOST CENTRALLY, i.e. with
the smallest `max|natural coordinate|`. On a genuine shared face the candidates agree
(C0), so this only matters at edges/overlaps, where it makes the result single-valued
and independent of bucket order (an all-interior neighbour no longer wins over the
surface-straddling element that truly owns the point).
"""
function eval_sdf(s::UnstructuredSDF, p::SVector{3,Float64})
    cell = hash_cell(s.bmin, s.inv_cell, s.nc, p)
    bucket = get(s.buckets, cell, nothing)
    if bucket !== nothing
        best_e = 0
        best_xi = zero(SVector{3,Float64})
        best_depth = Inf                       # smallest max|xi| among containing elements
        for e in bucket
            X = hex_nodes(s, e)
            xi, inside = inverse_hex_map(X, p)
            if inside
                depth = maximum(abs.(xi))
                if depth < best_depth
                    best_depth = depth
                    best_xi = xi
                    best_e = e
                end
            end
        end
        if best_e != 0
            N = hex8_shape(best_xi)
            hx = s.hexes[best_e]
            val = 0.0
            for a = 1:8
                val += N[a] * s.phi[hx[a]]
            end
            return val
        end
    end
    # Outside the meshed domain: positive distance to the domain box (see docstring).
    q = clamp.(p, s.bmin, s.bmax)
    return sqrt(sum(abs2, p - q))
end
