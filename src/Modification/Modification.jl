module Modification

export BoundedPlane, warp_mesh_by_planes_sdf!, PlaneDefinition

using StaticArrays
using LinearAlgebra

using Implicit2TetMesh.Fundamentals

include("CuttingPlaneTypes.jl")
include("ModifyResultingMesh.jl")

end
