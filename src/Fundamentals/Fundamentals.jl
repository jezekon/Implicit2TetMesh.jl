module Fundamentals

export BlockMesh, get_cell_sdf_values, eval_sdf, compute_gradient, calculate_volume_from_sdf

using StaticArrays
using LinearAlgebra
using FastGaussQuadrature
using Base.Threads

include("BlockMesh.jl")
include("SDFOperations.jl")
include("CalcVolumeFromSDF.jl")

end
