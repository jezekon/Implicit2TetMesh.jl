module Modification

export BoundedPlane,
    warp_mesh_by_planes_sdf!,
    PlaneDefinition,
    Rectangle,
    Square,
    Circle,
    Ellipse,
    remove_isolated_components!

using StaticArrays
using LinearAlgebra
using Printf

using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh

include("CuttingPlaneTypes.jl")
include("ModifyResultingMesh.jl")
include("RemoveIsolatedComponents.jl")

end
