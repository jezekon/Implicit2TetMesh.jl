# ----------------------------------------------------------------------------
# SDF source abstraction
# ----------------------------------------------------------------------------
# The mesher consumes its input field through exactly one scalar query, eval_sdf,
# plus a bounding box. Turning that query into an interface lets the SAME pipeline
# (lattice fill, warp, stencils, connectivity) mesh either a structured grid
# (trilinear interpolation -- the reference path) or an unstructured conforming
# HEX8 finite-element field (isoparametric shape functions), without any change to
# the meshing code. The generation lattice itself always stays structured -- that
# is a property of isosurface stuffing; only the SDF SOURCE is pluggable.
#
# CONVENTION (Etapa 1): phi < 0 = inside, phi > 0 = outside, phi = 0 = on the
# surface. Every concrete source must already store the field in this sign; the
# adapters in Adapters.jl convert foreign conventions (e.g. SIMP densities, or a
# level set with positive = inside) on the way in.

"""
    SDFSource

Abstract supertype for a pluggable implicit field. A concrete source implements
two methods:

  - `eval_sdf(src, p::SVector{3,Float64})::Float64` -- the field value at point
    `p` (phi < 0 inside, phi > 0 outside, phi = 0 on the surface).
  - `bbox(src) -> (min_corner, max_corner)` -- the axis-aligned bounding box of
    the region where the field is defined, as two `SVector{3,Float64}`. The mesher
    builds its lattice from this box plus a padding ring.

The field does NOT have to be a true signed distance function. After the
edge-based warp (Etapa 2) the cut point on an edge is `alpha = phi_i/(phi_i-phi_j)`
(scale-invariant) and no gradients are used, so only the zero set and the sign
near it matter -- SIMP density fields and level sets qualify once mapped to the
phi<0=inside convention.

Concrete sources: [`StructuredSDF`](@ref) (uniform grid, trilinear) and
[`UnstructuredSDF`](@ref) (conforming HEX8 mesh, FE shape functions).
"""
abstract type SDFSource end

"""
    bbox(src::SDFSource) -> (min_corner, max_corner)

Axis-aligned bounding box of the source's domain, returned as two
`SVector{3,Float64}` corners `(min, max)`. Every concrete `SDFSource` overrides
this method.
"""
function bbox end
