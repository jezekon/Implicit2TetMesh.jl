module GenerateMesh

export generate_mesh!, warp!, update_connectivity!, slice_ambiguous_tetrahedra!, adjust_nodes_to_isosurface!, export_mesh_vtk, warp_node_to_isocontour!, longest_edge

using Statistics
using StaticArrays
using LinearAlgebra
using WriteVTK

using Implicit2TetMesh.Fundamentals

include("Schemes/A15Scheme.jl")
include("Schemes/SchlafliScheme.jl")
include("TetGenerator.jl")
include("Stencils.jl")
include("ExportMesh.jl")


end
