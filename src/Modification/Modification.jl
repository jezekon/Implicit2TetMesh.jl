module Modification

export BoundedPlane,
    warp_mesh_by_planes_sdf!,
    PlaneDefinition,
    Rectangle,
    Square,
    Circle,
    Ellipse,
    remove_isolated_components!,
    RelaxOptions,
    relax_mesh!,
    CapRecoveryOptions,
    recover_boundary_caps!

using StaticArrays
using LinearAlgebra
using Printf
using Statistics

using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh
# The single EXACT orientation test lives in GenerateMesh and is intentionally not
# exported; the relaxation gate routes every orientation decision through it.
using Implicit2TetMesh.GenerateMesh: is_positively_oriented

include("CuttingPlaneTypes.jl")
include("ModifyResultingMesh.jl")
include("RemoveIsolatedComponents.jl")
include("RelaxMesh.jl")
# CapRecovery reuses RelaxMesh's face_outward_normal + reproject_to_surface and
# RemoveIsolatedComponents' face_key, so it is included after both.
include("CapRecovery.jl")

end
