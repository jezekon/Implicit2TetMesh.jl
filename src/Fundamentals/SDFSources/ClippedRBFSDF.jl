# ----------------------------------------------------------------------------
# RBF interpolant clipped by the design domain (domain clip of rho2sdf.jl)
# ----------------------------------------------------------------------------

"""
    ClippedRBFSDF(rbf::RBFSDF, domain::DomainDistance)

The smoothed SDF of rho2sdf.jl clipped by the design domain D. In the rho2sdf convention
(positive inside) the body is

    {min(f(x) - clip_shift, d_D(x)) >= 0} = {f >= clip_shift} ∩ D,

f the RBF sum of `rbf` (whose `shift` is clip_shift) and d_D the signed distance of D
(`DomainDistance`, positive inside). rho2sdf.jl sets clip_shift so that this clipped body keeps
the volume of the SIMP result. In this package's convention (phi < 0 inside) the source is

    phi(x) = max(eval_sdf(rbf, x), -d_D(x)).

`max` is the exact intersection of the two bodies; the field has a crease where the RBF surface
meets the boundary of D, so mesh it with `cut_points = :bisection`, which places the surface
vertices on the zero set itself. Outside the RBF grid box the source returns the distance to the
box, as `RBFSDF`. `h` is the RBF grid step (`dx = source.h`, as for `RBFSDF`).

Built by [`load_rbf_source`](@ref) from a JLD2 written with `domain_clip = true`.
"""
struct ClippedRBFSDF <: SDFSource
    rbf::RBFSDF                   # with shift = clip_shift
    domain::DomainDistance        # d_D, positive inside
    h::Float64                    # RBF grid step (= rbf.h)
end

ClippedRBFSDF(rbf::RBFSDF, domain::DomainDistance) = ClippedRBFSDF(rbf, domain, rbf.h)

bbox(s::ClippedRBFSDF) = bbox(s.rbf)

"""
    eval_sdf(s::ClippedRBFSDF, p) -> Float64

max(RBF field, -d_D) at `p` (phi < 0 inside); outside the grid box the distance to the box.
"""
function eval_sdf(s::ClippedRBFSDF, p::SVector{3,Float64})
    if any(p .< s.rbf.bmin) || any(p .> s.rbf.bmax)
        return eval_sdf(s.rbf, p)
    end
    return max(eval_sdf(s.rbf, p), -s.domain(p))
end
