# ----------------------------
# Helper function: Get SDF values at the 8 corners of a cell (unchanged)
# ----------------------------
function get_cell_sdf_values(mesh::BlockMesh, i::Int, j::Int, k::Int)
    1 <= i < mesh.nx || throw(BoundsError(mesh.SDF, i))
    1 <= j < mesh.ny || throw(BoundsError(mesh.SDF, j))
    1 <= k < mesh.nz || throw(BoundsError(mesh.SDF, k))
    return SVector{8}(
        mesh.SDF[i, j, k],       # front-bottom-left
        mesh.SDF[i+1, j, k],     # front-bottom-right
        mesh.SDF[i+1, j+1, k],   # front-top-right
        mesh.SDF[i, j+1, k],     # front-top-left
        mesh.SDF[i, j, k+1],     # back-bottom-left
        mesh.SDF[i+1, j, k+1],   # back-bottom-right
        mesh.SDF[i+1, j+1, k+1], # back-top-right
        mesh.SDF[i, j+1, k+1],    # back-top-left
    )
end


# ----------------------------
# Helper function: SDF evaluation, delegated to the pluggable field source
# ----------------------------
# The whole pipeline (warp node_sdf, bisection cut point, surface-tet centroid test)
# queries the field ONLY through this function, so making it a thin delegation to
# mesh.sdf_source is the single seam that lets the mesher consume either a structured
# grid (trilinear) or an unstructured HEX8 field (FE shape functions). For structured
# input the source is a StructuredSDF wrapping the same grid + values, so this returns
# exactly the original trilinear interpolation (bit-identical). The lattice-corner
# cache (get_cell_sdf_values above) keeps reading mesh.SDF directly.
function eval_sdf(mesh::BlockMesh, p::SVector{3,Float64})
    return eval_sdf(mesh.sdf_source, p)
end

