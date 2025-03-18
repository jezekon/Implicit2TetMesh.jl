module Modification

export BoundedPlane, warp_mesh_by_planes_sdf!, PlaneDefinition, Rectangle, Square, Circle, Ellipse

using StaticArrays
using LinearAlgebra

using Implicit2TetMesh.Fundamentals

include("CuttingPlaneTypes.jl")
include("ModifyResultingMesh.jl")

end
