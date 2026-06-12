# ----------------------------------------------------------------------------
# Structured-grid SDF source (trilinear interpolation)
# ----------------------------------------------------------------------------

"""
    StructuredSDF(grid, values)

Implicit field sampled on a uniform structured grid; between grid points the field
is the trilinear interpolation of the eight surrounding corner values. This is the
reference path -- it reproduces the original `eval_sdf` exactly -- and the fast
path the auto-detector ([`build_sdf_source`](@ref)) falls back to for
tensor-product inputs.

# Arguments
- `grid::Array{SVector{3,Float64},3}`: node coordinates of the grid.
- `values::Array{Float64,3}`: field at those nodes, already in the phi<0=inside
  convention. Must have the same size as `grid`.

The grid must be uniform along each axis: a query point is mapped to fractional
indices through the box extents `grid[1,1,1]` / `grid[end,end,end]`, exactly as the
original structured `eval_sdf` did.
"""
struct StructuredSDF <: SDFSource
    vmin::SVector{3,Float64}      # grid[1, 1, 1]
    vmax::SVector{3,Float64}      # grid[end, end, end]
    nx::Int
    ny::Int
    nz::Int
    values::Array{Float64,3}      # corner field values, phi < 0 = inside
end

function StructuredSDF(grid::Array{SVector{3,Float64},3}, values::Array{Float64,3})
    size(grid) == size(values) || error(
        "StructuredSDF: grid size $(size(grid)) does not match values size $(size(values))",
    )
    nx, ny, nz = size(grid)
    return StructuredSDF(grid[1, 1, 1], grid[end, end, end], nx, ny, nz, values)
end

"""
    bbox(s::StructuredSDF) -> (min_corner, max_corner)

Bounding box of the structured grid (its first and last node).
"""
bbox(s::StructuredSDF) = (s.vmin, s.vmax)

"""
    eval_sdf(s::StructuredSDF, p) -> Float64

Trilinear interpolation of the grid field at `p`. The query point is mapped to
fractional grid indices through the box extents and clamped to the last cell, so
points on or just outside the grid border are evaluated against the boundary cell.
This is the original structured `eval_sdf` body, moved here verbatim so the
structured path stays bit-identical.
"""
function eval_sdf(s::StructuredSDF, p::SVector{3,Float64})
    # Minimum and maximum grid coordinates.
    vmin = s.vmin
    vmax = s.vmax

    # Normalize point coordinates to the interval [0, 1].
    r = (p .- vmin) ./ (vmax .- vmin)

    # Convert to (fractional) grid indices.
    i_f = r[1] * (s.nx - 1) + 1
    j_f = r[2] * (s.ny - 1) + 1
    k_f = r[3] * (s.nz - 1) + 1

    i0 = clamp(floor(Int, i_f), 1, s.nx - 1)
    j0 = clamp(floor(Int, j_f), 1, s.ny - 1)
    k0 = clamp(floor(Int, k_f), 1, s.nz - 1)

    i1 = i0 + 1
    j1 = j0 + 1
    k1 = k0 + 1

    # Local weights.
    xd = i_f - i0
    yd = j_f - j0
    zd = k_f - k0

    # SDF values at the eight corners of the cell.
    c000 = s.values[i0, j0, k0]
    c100 = s.values[i1, j0, k0]
    c010 = s.values[i0, j1, k0]
    c110 = s.values[i1, j1, k0]
    c001 = s.values[i0, j0, k1]
    c101 = s.values[i1, j0, k1]
    c011 = s.values[i0, j1, k1]
    c111 = s.values[i1, j1, k1]

    # Trilinear interpolation.
    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    return c0 * (1 - zd) + c1 * zd
end
