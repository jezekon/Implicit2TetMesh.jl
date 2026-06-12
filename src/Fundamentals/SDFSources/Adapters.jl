# ----------------------------------------------------------------------------
# Input adapters: build an SDFSource from topology-optimization / level-set data
# ----------------------------------------------------------------------------
# These turn the two common unstructured inputs -- SIMP element densities and a
# nodal level set -- into the phi<0=inside field the mesher consumes, and pick the
# concrete source type (structured fast path vs unstructured HEX8) automatically.

# Tolerance for recognizing coincident coordinates when auto-detecting a grid;
# derived from the data so it scales with the model size.
function default_grid_tol(nodes::Vector{SVector{3,Float64}})
    bmin = nodes[1]
    bmax = nodes[1]
    for p in nodes
        bmin = min.(bmin, p)
        bmax = max.(bmax, p)
    end
    span = maximum(bmax - bmin)
    return max(1e-12, 1e-6 * span)
end

# Sorted distinct values of a coordinate, merging entries closer than `tol`.
function unique_sorted_coords(vals::Vector{Float64}, tol::Float64)
    s = sort(vals)
    u = Float64[s[1]]
    for v in s
        if v - u[end] > tol
            push!(u, v)
        end
    end
    return u
end

# True if the sorted coordinate levels are equally spaced (within `tol`).
function is_uniform_spacing(u::Vector{Float64}, tol::Float64)
    length(u) <= 2 && return true
    d0 = u[2] - u[1]
    for i = 2:(length(u)-1)
        abs((u[i+1] - u[i]) - d0) > tol && return false
    end
    return true
end

# Index of the level matching `x` within `tol`, or 0 if none.
function coord_index(u::Vector{Float64}, x::Float64, tol::Float64)
    i = searchsortedfirst(u, x - tol)
    (1 <= i <= length(u) && abs(u[i] - x) <= tol) && return i
    return 0
end

"""
    detect_structured(nodes, phi, tol) -> Union{StructuredSDF, Nothing}

Recognize a UNIFORM tensor-product grid from a node set: the nodes must be the full
Cartesian product of equally-spaced x/y/z levels (within `tol`). On success returns
the equivalent [`StructuredSDF`](@ref) (the fast path); otherwise `nothing` (a graded
or genuinely unstructured mesh, handled by [`UnstructuredSDF`](@ref)).
"""
function detect_structured(nodes::Vector{SVector{3,Float64}}, phi::Vector{Float64}, tol::Float64)
    ux = unique_sorted_coords([p[1] for p in nodes], tol)
    uy = unique_sorted_coords([p[2] for p in nodes], tol)
    uz = unique_sorted_coords([p[3] for p in nodes], tol)
    nx, ny, nz = length(ux), length(uy), length(uz)

    # Must be a complete, uniformly spaced lattice.
    nx * ny * nz == length(nodes) || return nothing
    (is_uniform_spacing(ux, tol) && is_uniform_spacing(uy, tol) && is_uniform_spacing(uz, tol)) ||
        return nothing

    grid = Array{SVector{3,Float64},3}(undef, nx, ny, nz)
    values = Array{Float64,3}(undef, nx, ny, nz)
    filled = falses(nx, ny, nz)
    for (n, p) in enumerate(nodes)
        ix = coord_index(ux, p[1], tol)
        iy = coord_index(uy, p[2], tol)
        iz = coord_index(uz, p[3], tol)
        (ix == 0 || iy == 0 || iz == 0) && return nothing
        filled[ix, iy, iz] && return nothing        # duplicate node -> not a clean grid
        filled[ix, iy, iz] = true
        grid[ix, iy, iz] = p
        values[ix, iy, iz] = phi[n]
    end
    all(filled) || return nothing
    return StructuredSDF(grid, values)
end

"""
    build_sdf_source(nodes, hexes, nodal_phi; force = :auto, grid_tol = nothing) -> SDFSource

Pick the concrete source for a HEX8 field given node coordinates, element
connectivity, and a nodal field (already phi<0=inside):

  - `:auto` (default): a uniform tensor-product grid becomes a [`StructuredSDF`](@ref)
    (fast trilinear path); anything else becomes an [`UnstructuredSDF`](@ref).
  - `:structured`: force the grid fast path (errors if the nodes are not a uniform grid).
  - `:unstructured`: force the HEX8 FE path (used by the round-trip test to exercise
    point location + inverse mapping even on grid-derived inputs).

`grid_tol` overrides the coordinate-merging tolerance of the auto-detector.
"""
function build_sdf_source(nodes, hexes, nodal_phi; force::Symbol = :auto, grid_tol = nothing)
    nodes_sv = SVector{3,Float64}[SVector{3,Float64}(p[1], p[2], p[3]) for p in nodes]
    phi_f = Float64.(nodal_phi)
    tol = grid_tol === nothing ? default_grid_tol(nodes_sv) : Float64(grid_tol)

    if force === :unstructured
        return UnstructuredSDF(nodes_sv, hexes, phi_f)
    elseif force === :structured
        s = detect_structured(nodes_sv, phi_f, tol)
        s === nothing &&
            error("build_sdf_source(force = :structured): nodes are not a uniform tensor-product grid")
        return s
    elseif force === :auto
        s = detect_structured(nodes_sv, phi_f, tol)
        return s === nothing ? UnstructuredSDF(nodes_sv, hexes, phi_f) : s
    else
        error("build_sdf_source: force must be :auto, :structured or :unstructured (got $force)")
    end
end

"""
    simp_to_sdf_source(nodes, hexes, densities; iso_level = 0.5, force = :auto,
                       grid_tol = nothing) -> SDFSource

Adapter for a SIMP topology-optimization result. `densities` are element-constant
(one per hex). They are averaged to the nodes (volume-weighted over the incident
elements) and turned into the field `phi = iso_level - rho`, so `phi < 0` exactly
where the smoothed density exceeds `iso_level` (solid). `iso_level` defaults to the
usual 0.5 density threshold and is user-adjustable. The nodal field is then handed
to [`build_sdf_source`](@ref).
"""
function simp_to_sdf_source(
    nodes,
    hexes,
    densities;
    iso_level::Real = 0.5,
    force::Symbol = :auto,
    grid_tol = nothing,
)
    nodes_sv = SVector{3,Float64}[SVector{3,Float64}(p[1], p[2], p[3]) for p in nodes]
    length(densities) == length(hexes) ||
        error("simp_to_sdf_source: $(length(densities)) densities for $(length(hexes)) elements")

    nnodes = length(nodes_sv)
    weighted = zeros(Float64, nnodes)     # sum of vol_e * rho_e over incident elements
    weight = zeros(Float64, nnodes)       # sum of vol_e
    for (e, h) in enumerate(hexes)
        X = ntuple(a -> nodes_sv[Int(h[a])], 8)
        ve = hex_volume(X)
        re = Float64(densities[e])
        for a = 1:8
            nd = Int(h[a])
            weighted[nd] += ve * re
            weight[nd] += ve
        end
    end

    phi = Vector{Float64}(undef, nnodes)
    for i = 1:nnodes
        rho = weight[i] > 0 ? weighted[i] / weight[i] : 0.0
        phi[i] = iso_level - rho
    end
    return build_sdf_source(nodes_sv, hexes, phi; force = force, grid_tol = grid_tol)
end

"""
    levelset_to_sdf_source(nodes, hexes, nodal_values; inside = :negative,
                           force = :auto, grid_tol = nothing) -> SDFSource

Adapter for a nodal level-set field. `inside` states the input sign convention so
the stored field ends up phi<0=inside (Etapa 1): `:negative` (the field is already
negative inside -- used as is) or `:positive` (positive inside -- negated). The
field is then handed to [`build_sdf_source`](@ref).
"""
function levelset_to_sdf_source(
    nodes,
    hexes,
    nodal_values;
    inside::Symbol = :negative,
    force::Symbol = :auto,
    grid_tol = nothing,
)
    phi = if inside === :negative
        Float64.(nodal_values)
    elseif inside === :positive
        -Float64.(nodal_values)
    else
        error("levelset_to_sdf_source: inside must be :negative or :positive (got $inside)")
    end
    return build_sdf_source(nodes, hexes, phi; force = force, grid_tol = grid_tol)
end
