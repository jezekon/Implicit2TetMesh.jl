"""
    calculate_volume_from_sdf(
        fine_sdf::Array{Float32,3}, 
        fine_grid::Array{Vector{Float32},3}; 
        iso_threshold::Float32=0.0f0, 
        detailed_quad_order::Int=9
    ) -> Float32

Calculates volume of geometry defined by SDF iso-surface.

# Arguments
- `fine_sdf`: SDF values at regular grid nodes
- `fine_grid`: Coordinates of grid nodes, each element is [x,y,z]
- `iso_threshold`: Iso-surface value (default = 0.0)
- `detailed_quad_order`: Gauss-Legendre quadrature order

# Returns
- Volume of the geometry

# Description
Assumes the phi < 0 = inside convention. For each grid element:
1. Skip if outside iso-surface (all SDF values > threshold)
2. Add full volume if inside iso-surface (all SDF values ≤ threshold)
3. Use quadrature for partial elements intersecting the iso-surface
"""
function calculate_volume_from_sdf(
    fine_sdf::Array{Float32,3},
    fine_grid::Array{Vector{Float32},3};
    iso_threshold::Float32 = 0.0f0,
    detailed_quad_order::Int = 9,
)
    # Determine grid dimensions
    nx, ny, nz = size(fine_sdf)
    @assert size(fine_grid) == (nx, ny, nz) "Dimensions of fine_sdf and fine_grid must match"

    # Calculate element size (assuming all elements are cubes with the same edge length)
    edge_vector = fine_grid[2, 1, 1] .- fine_grid[1, 1, 1]
    element_edge_length = norm(edge_vector)
    element_volume = element_edge_length^3

    # Set up Gauss-Legendre quadrature
    gp_double, w_double = FastGaussQuadrature.gausslegendre(detailed_quad_order)
    gp = Float32.(gp_double)
    w = Float32.(w_double)

    # Initialize total volume using atomic variable for parallel processing
    total_volume = Atomic{Float32}(0.0f0)

    # Precompute Jacobian for cubic element
    # For mapping from [-1,1]^3 to a cube with edge a, the Jacobian determinant is (a/2)^3 = a^3/8
    jacobian_det = element_volume / 8.0f0

    # Parallel iteration over all elements
    @threads for k = 1:(nz-1)
        for j = 1:(ny-1)
            for i = 1:(nx-1)
                # Get SDF values at element nodes
                c000 = fine_sdf[i, j, k]
                c100 = fine_sdf[i+1, j, k]
                c010 = fine_sdf[i, j+1, k]
                c110 = fine_sdf[i+1, j+1, k]
                c001 = fine_sdf[i, j, k+1]
                c101 = fine_sdf[i+1, j, k+1]
                c011 = fine_sdf[i, j+1, k+1]
                c111 = fine_sdf[i+1, j+1, k+1]

                # Check if element can contain the iso-surface
                min_value = min(c000, c100, c010, c110, c001, c101, c011, c111)
                max_value = max(c000, c100, c010, c110, c001, c101, c011, c111)

                # Skip element if completely outside the iso-surface
                if min_value > iso_threshold
                    continue
                end

                # If element is completely inside the iso-surface, add its full volume
                if max_value <= iso_threshold
                    atomic_add!(total_volume, element_volume)
                    continue
                end

                # Element intersects the iso-surface, perform numerical integration
                element_volume_partial = 0.0f0

                # Perform numerical integration over all quadrature points
                for kq = 1:detailed_quad_order
                    zeta_m = gp[kq]  # Gauss point in [-1,1]
                    zeta = (zeta_m + 1) / 2  # Convert to [0,1] for interpolation

                    for jq = 1:detailed_quad_order
                        eta_m = gp[jq]
                        eta = (eta_m + 1) / 2

                        for iq = 1:detailed_quad_order
                            xi_m = gp[iq]
                            xi = (xi_m + 1) / 2

                            # Trilinear interpolation of SDF value at Gauss point
                            c00 = c000 * (1.0f0 - xi) + c100 * xi
                            c01 = c001 * (1.0f0 - xi) + c101 * xi
                            c10 = c010 * (1.0f0 - xi) + c110 * xi
                            c11 = c011 * (1.0f0 - xi) + c111 * xi

                            c0 = c00 * (1.0f0 - eta) + c10 * eta
                            c1 = c01 * (1.0f0 - eta) + c11 * eta

                            point_sdf = c0 * (1.0f0 - zeta) + c1 * zeta

                            # Add contribution only if point is inside (or on) the iso-surface
                            if point_sdf <= iso_threshold
                                weight = w[iq] * w[jq] * w[kq]
                                element_volume_partial += weight * jacobian_det
                            end
                        end
                    end
                end

                # Add to total volume
                atomic_add!(total_volume, element_volume_partial)
            end
        end
    end

    return total_volume[]
end

"""
    calculate_volume_from_sdf(source::UnstructuredSDF; iso_threshold=0.0,
                              quad_order=9) -> Float64

Volume of the solid region (`phi <= iso_threshold`) of an unstructured HEX8 field --
the per-source counterpart of the structured array method above (Etapa 8). The field
inside each hex is the trilinear FE shape-function interpolation, which is a convex
combination of the 8 nodal values, so an element is fully solid when all its nodal
values are `<= iso_threshold` and fully void when all are `> iso_threshold`. Those two
cases are handled in closed form (`hex_volume` is exact for the trilinear map); only
boundary elements are integrated:

  - fully solid  -> add `hex_volume(X)`;
  - fully void   -> skip;
  - boundary     -> `quad_order^3` Gauss-Legendre points; add `w * |det J|` at each
    point whose interpolated field is `<= iso_threshold`.

`quad_order` controls only how finely the cut surface is resolved inside a boundary
element (the integrand is a step function there); it defaults to 9 to match the
structured method.
"""
function calculate_volume_from_sdf(
    source::UnstructuredSDF;
    iso_threshold::Float64 = 0.0,
    quad_order::Int = 9,
)
    quad_order >= 1 || error("calculate_volume_from_sdf: quad_order must be >= 1, got $quad_order")
    gp, w = FastGaussQuadrature.gausslegendre(quad_order)

    total_volume = Atomic{Float64}(0.0)
    @threads for e = 1:length(source.hexes)
        hx = source.hexes[e]
        X = ntuple(a -> source.nodes[hx[a]], 8)
        pv = ntuple(a -> source.phi[hx[a]], 8)
        mn, mx = extrema(pv)

        # Fully void: every interior value exceeds the threshold.
        mn > iso_threshold && continue

        # Fully solid: the whole element is below the threshold -> exact hex volume.
        if mx <= iso_threshold
            atomic_add!(total_volume, hex_volume(X))
            continue
        end

        # Boundary element: integrate the indicator of {phi <= iso} times |det J|.
        partial = 0.0
        for kq = 1:quad_order, jq = 1:quad_order, iq = 1:quad_order
            xi = SVector(gp[iq], gp[jq], gp[kq])
            N = hex8_shape(xi)
            phi_q = 0.0
            for a = 1:8
                phi_q += N[a] * pv[a]
            end
            if phi_q <= iso_threshold
                _, J = hex8_map_and_jacobian(X, xi)
                partial += w[iq] * w[jq] * w[kq] * abs(det(J))
            end
        end
        atomic_add!(total_volume, partial)
    end
    return total_volume[]
end
