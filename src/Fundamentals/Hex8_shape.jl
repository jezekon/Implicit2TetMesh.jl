# ----------------------------
# Helper function to compute shape functions for a hex8 element
# Accepts local coordinates in [0,1] (normalized tile_ref) and returns shape functions
# computed with standard transformation to [-1,1]
# ----------------------------
function shape_functions(ξηζ::SVector{3,Float64})::SVector{8,Float64}
    # Transform local coordinates from [0,1] to [-1,1]
    ξ = 2 * ξηζ[1] - 1.0
    η = 2 * ξηζ[2] - 1.0
    ζ = 2 * ξηζ[3] - 1.0
    coef = 1 / 8.0
    return @SVector [
        coef * (1 - ξ) * (1 - η) * (1 - ζ),
        coef * (1 + ξ) * (1 - η) * (1 - ζ),
        coef * (1 + ξ) * (1 + η) * (1 - ζ),
        coef * (1 - ξ) * (1 + η) * (1 - ζ),
        coef * (1 - ξ) * (1 - η) * (1 + ζ),
        coef * (1 + ξ) * (1 - η) * (1 + ζ),
        coef * (1 + ξ) * (1 + η) * (1 + ζ),
        coef * (1 - ξ) * (1 + η) * (1 + ζ),
    ]
end
