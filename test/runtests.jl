using Test
using JLD2
using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh
using Implicit2TetMesh.Optimization
using Implicit2TetMesh.Modification
using Implicit2TetMesh.Utils


@testset "Implicit2TetMesh.jl" begin

  RUN_beam = true
  RUN_main = false
  RUN_main_param = false

  if RUN_beam

      # taskName = "beam_vfrac_04_B-1.0_smooth-1_Approx"
      taskName = "beam_vfrac_04_B-1.0_smooth-1_Interp"
      # taskName = "beam_vfrac_04_B-1.0_smooth-2_Approx"
      # taskName = "beam_vfrac_04_B-1.0_smooth-2_Interp"

      # @load "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Approximation.jld2" fine_grid
      # @load "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Approximation.jld2" fine_sdf
      @load "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2" fine_grid
      @load "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2" fine_sdf
      # @load "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-2_Approximation.jld2" fine_grid
      # @load "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-2_Approximation.jld2" fine_sdf
      # @load "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-2_Interpolation.jld2" fine_grid
      # @load "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-2_Interpolation.jld2" fine_sdf
 
      scheme = "A15"
      # scheme = "Schlafli"
      
      # Definice rovin
      plane_definitions = [
        PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.)),
        PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.))
      ]
      warp_param = 0.3

      mesh = BlockMesh(fine_sdf, fine_grid)

      # Choose scheme: "A15" or "Schlafli"
      generate_mesh!(mesh, scheme)

      warp!(mesh, scheme) # Warp nearest nodes to the isocontour
      # export_mesh_vtu_quality(mesh, "$(taskName)_TriMesh-step_warp$(scheme).vtu")

      slice_ambiguous_tetrahedra!(mesh) # Remove elements outside the body
      # adjust_nodes_to_isosurface!(mesh) # Simple cut of elements to follow the isocontour
 
      update_connectivity!(mesh)

      export_mesh_vtu_quality(mesh, "$(taskName)_TriMesh-notOPT_$(scheme).vtu")

      TetMesh_volumes(mesh)

      remove_inverted_elements!(mesh)
      update_connectivity!(mesh)
      optimize_mesh!(mesh, scheme)

      # Korekce objemu
      success = correct_mesh_volume!(mesh, fine_sdf, fine_grid)
      
      if success
          @info "Volume correction successful"
      else
          @warn "Volume correction did not converge to required tolerance"
      end
      
      # Vyhodnocení přesnosti
      assess_volume_accuracy(mesh, fine_sdf, fine_grid)

      export_mesh_vtu(mesh, "$(taskName)_TriMesh-$(scheme).vtu")

      # Apply cutting planes only if they are defined
      if !isempty(plane_definitions)
        warp_mesh_by_planes_sdf!(mesh, plane_definitions, warp_param)
        update_connectivity!(mesh)
        export_mesh_vtu(mesh, "$(taskName)_TriMesh-$(scheme)_cut.vtu")
      end

      TetMesh_volumes(mesh)
 
      slice_mesh_with_plane!(mesh, "x", 0.977, -1)
      export_mesh_vtu(mesh, "$(taskName)_TriMesh-$(scheme)_plane.vtu")
      # assess_mesh_quality(mesh, "mesh_quality")
  end

  if RUN_main
    mesh = generate_tetrahedral_mesh(
    "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2",
    "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2",
      "cantilever_beam_interp")
  end

  if RUN_main_param

    plane_definitions = [
        PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
        PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0))]
    
    options=MeshGenerationOptions(
            scheme = "A15",
            optimize = true,
            split_elements = false,
            quality_export = true,
            warp_param = 0.,
            plane_definitions = plane_definitions)

    mesh = generate_tetrahedral_mesh(
        "../data/Z_cantilever_beam_vfrac_04_FineGrid_B-1.0_smooth-1_Interpolation.jld2",
        "../data/Z_cantilever_beam_vfrac_04_FineSDF_B-1.0_smooth-1_Interpolation.jld2",
        "cantilever_beam_interp_cut",
        options=options)

    # assess_mesh_quality(mesh, "mesh_quality")
  end
end
