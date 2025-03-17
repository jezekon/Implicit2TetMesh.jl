module Fundamentals

export BlockMesh, get_cell_sdf_values, eval_sdf, compute_gradient

using StaticArrays
using LinearAlgebra

include("BlockMesh.jl")
include("SDFOperations.jl")

end
