# ----------------------------------------------------------------------------
# Gaussian RBF interpolant as an SDF source (evaluated directly, no resampling)
# ----------------------------------------------------------------------------

"""
    RBFSDF(weights, origin, h; sigma = h, threshold = 1e-3, shift = 0.0, inside = :positive)

Smoothed SDF given as a Gaussian RBF interpolant with centres on a regular grid (nodes
`origin + h·(i, j, k)`, `i = 0 … size(weights, 1) - 1`, …):

    f(x) = Σ_i w_i exp(-(|x - x_i| / sigma)²) - shift,

with every kernel cut off where its value is ≤ `threshold` (the truncation used to assemble
the interpolation matrix in rho2sdf.jl). `eval_sdf` evaluates this sum at the query point
itself, so the mesher sees the C∞ field at its lattice nodes and, with
`cut_points = :bisection`, at the true zeros along lattice edges. No trilinear resampling
is involved, and the lattice can keep the spacing `h` of the RBF grid.

`inside = :positive` (rho2sdf.jl convention) negates the field into this package's
phi < 0 = inside convention; `:negative` uses it as is. Outside the grid box the source
returns the (positive) distance to the box.
"""
struct RBFSDF <: SDFSource
    weights::Array{Float64,3}
    origin::SVector{3,Float64}
    h::Float64
    a::Float64                    # exponent factor in grid units, (h / sigma)²
    threshold::Float64
    reach::Float64                # cut-off radius in grid units
    shift::Float64
    sgn::Float64                  # +1 or -1: maps the input convention to phi < 0 = inside
    bmin::SVector{3,Float64}
    bmax::SVector{3,Float64}
end

function RBFSDF(
    weights::AbstractArray{<:Real,3},
    origin::AbstractVector{<:Real},
    h::Real;
    sigma::Real = h,
    threshold::Real = 1e-3,
    shift::Real = 0.0,
    inside::Symbol = :positive,
)
    inside in (:positive, :negative) ||
        error("RBFSDF: inside must be :positive or :negative, got :$inside")
    0 < threshold < 1 || error("RBFSDF: threshold must lie in (0, 1), got $threshold")
    o = SVector{3,Float64}(origin...)
    n = size(weights)
    bmax = o .+ h .* SVector{3,Float64}(n[1] - 1, n[2] - 1, n[3] - 1)
    return RBFSDF(
        Float64.(weights), o, Float64(h), (h / sigma)^2, Float64(threshold),
        (sigma / h) * sqrt(-log(threshold)), Float64(shift),
        inside === :positive ? -1.0 : 1.0, o, bmax,
    )
end

bbox(s::RBFSDF) = (s.bmin, s.bmax)

# Index window of the kernels that reach the grid coordinate q along one axis
@inline function rbf_window(q::Float64, reach::Float64, n::Int)
    return max(0, ceil(Int, q - reach)), min(n - 1, floor(Int, q + reach))
end

"""
    eval_sdf(s::RBFSDF, p) -> Float64

Value of the RBF interpolant at `p` (sum over the kernels whose cut-off reaches `p`),
shifted and mapped to phi < 0 = inside. Outside the grid box: distance to the box.
"""
function eval_sdf(s::RBFSDF, p::SVector{3,Float64})
    if any(p .< s.bmin) || any(p .> s.bmax)
        q = clamp.(p, s.bmin, s.bmax)
        return sqrt(sum(abs2, p - q))
    end
    nx, ny, nz = size(s.weights)
    q = (p - s.origin) / s.h
    i0, i1 = rbf_window(q[1], s.reach, nx)
    j0, j1 = rbf_window(q[2], s.reach, ny)
    k0, k1 = rbf_window(q[3], s.reach, nz)

    # Separable kernel: exp(-a r²) = ex(i) ey(j) ez(k); the spherical cut-off is the test g > threshold
    ex = MVector{8,Float64}(undef)
    for i = i0:i1
        ex[i-i0+1] = exp(-s.a * (i - q[1])^2)
    end
    acc = 0.0
    @inbounds for k = k0:k1
        gz = exp(-s.a * (k - q[3])^2)
        for j = j0:j1
            gyz = gz * exp(-s.a * (j - q[2])^2)
            gyz > s.threshold || continue
            for i = i0:i1
                g = gyz * ex[i-i0+1]
                g > s.threshold && (acc += s.weights[i+1, j+1, k+1] * g)
            end
        end
    end
    return s.sgn * (acc - s.shift)
end
