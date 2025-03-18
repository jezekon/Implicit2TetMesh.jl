using Test
using JLD2
using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh
using Implicit2TetMesh.Optimization
using Implicit2TetMesh.Modification
using Implicit2TetMesh.Utils


@testset "Implicit2TetMesh.jl" begin

  # taskName = "beam_vfrac_04_B-1.0_smooth-1_Approximation"
  taskName = "beam_vfrac_04_B-1.0_smooth-1_Interpolation"

  # @load "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Approximation.jld2" fine_grid
  # @load "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Approximation.jld2" fine_sdf
  @load "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2" fine_grid
  @load "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2" fine_sdf
  
  scheme = "A15"
  
  # Definice rovin
  plane_definitions = [
    PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.)),
    PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.))
  ]
  warp_param = 0.3

  mesh = BlockMesh(fine_sdf, fine_grid)

  # Choose scheme: "A15" or "Schlafli"
  generate_mesh!(mesh, scheme)

  # assess_mesh_quality(mesh, "initial_mesh")
  warp!(mesh) # Warp nearest nodes to the isocontour

  update_connectivity!(mesh) # Update mesh topology

  slice_ambiguous_tetrahedra!(mesh) # Remove elements outside the body

  update_connectivity!(mesh)

  adjust_nodes_to_isosurface!(mesh) # Simple cut of elements to follow the isocontour

  TetMesh_volumes(mesh)
  optimize_mesh!(mesh)

  export_mesh_vtk(mesh, "$(taskName)_TriMesh.vtu")

  # Apply cutting planes only if they are defined
  if !isempty(plane_definitions)
    warp_mesh_by_planes_sdf!(mesh, plane_definitions, warp_param)
    update_connectivity!(mesh)
    export_mesh_vtk(mesh, "$(taskName)_TriMesh_cut.vtu")
  end

  TetMesh_volumes(mesh)
end
