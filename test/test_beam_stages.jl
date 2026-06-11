# ==============================================================================
# Per-stage beam diagnostic
# ==============================================================================
# Runs the beam pipeline stage by stage, mirroring the EXACT call sequence of
# generate_tetrahedral_mesh (src/TetMeshGenerator.jl), including the two
# update_connectivity!(build_ine = false) calls and the final build_ine = true
# refresh. After EVERY stage it FIRST exports the mesh to test/output/ with a
# numbered name, THEN runs that stage's assertions -- so when a stage breaks, the
# failing assertion names the phase and the VTU shows the damage.
#
# If the pipeline sequence in TetMeshGenerator.jl ever changes, mirror it here.
#
# Runnable standalone for interactive debugging:
#   julia --project=. test/test_beam_stages.jl
# or as part of the suite (included from runtests.jl).

using Test
using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh
using Implicit2TetMesh.Modification

# Load helpers + baselines only if not already present (idempotent include guard,
# so this works both standalone and when included after runtests.jl loaded them).
if !isdefined(@__MODULE__, :check_watertight)
    include(joinpath(@__DIR__, "helpers", "helpers.jl"))
end
if !isdefined(@__MODULE__, :BEAM_DIH)
    include(joinpath(@__DIR__, "helpers", "baselines.jl"))
end

@testset "Beam per-stage diagnostics" begin
    outdir = ensure_output_dir()
    fine_grid, fine_sdf = load_beam()
    mesh = BlockMesh(fine_sdf, fine_grid)

    # --- Stage 1: base A15 lattice -------------------------------------------
    @testset "Stage 1 generate_mesh!" begin
        generate_mesh!(mesh, "A15")
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_1-generated.vtu"))
        @test length(mesh.X) > 0
        @test length(mesh.IEN) > 0
        @test all_indices_in_bounds(mesh)
        @test all_finite_coords(mesh)
        @test count_inverted_exact(mesh) == 0
    end

    # --- Stage 2: edge-based warp onto the surface ---------------------------
    @testset "Stage 2 warp!" begin
        warp!(mesh, "A15")
        update_connectivity!(mesh; build_ine = false)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_2-warped.vtu"))
        @test all_finite_coords(mesh)
        nsc = node_sdf_consistency(mesh)
        @test nsc.pristine_within_tol     # untouched lattice nodes still match eval_sdf exactly
        @test nsc.n_on_surface > 0        # warp snapped some nodes onto the surface (node_sdf = 0)
    end

    # --- Stage 3: trim spikes / slice crossing tets --------------------------
    @testset "Stage 3 slice_ambiguous_tetrahedra!" begin
        slice_ambiguous_tetrahedra!(mesh, "A15")
        update_connectivity!(mesh; build_ine = false)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_3-sliced.vtu"))
        @test check_watertight(mesh).open_edges == 0
        @test all_indices_in_bounds(mesh)
    end

    # --- Stage 4: drop inverted elements -------------------------------------
    @testset "Stage 4 remove_inverted_elements!" begin
        remove_inverted_elements!(mesh)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_4-inverted_removed.vtu"))
        @test count_inverted_exact(mesh) == 0
    end

    # --- Stage 5: keep the largest component (final refresh builds INE) ------
    # This stage equals the no-planes pipeline output, so the dihedral baseline
    # is shared with the core suite (BEAM_DIH).
    @testset "Stage 5 remove_isolated_components!" begin
        remove_isolated_components!(mesh, keep_largest = true)
        update_connectivity!(mesh)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_5-components.vtu"))
        @test count_components(mesh) == 1
        @test check_watertight(mesh).open_edges == 0
        s = dihedral_stats(mesh)
        @test round(s.min; digits = 3) == BEAM_DIH.min
        @test round(s.max; digits = 3) == BEAM_DIH.max
        @test s.lt5 == BEAM_DIH.lt5
        @test s.lt10 == BEAM_DIH.lt10
        @test s.gt140 == BEAM_DIH.gt140
        @test s.inverted == BEAM_DIH.inverted
    end

    # --- Stage 6: warp surface nodes onto the cutting planes -----------------
    @testset "Stage 6 warp_mesh_by_planes_sdf!" begin
        plane_definitions = [
            PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
            PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0)),
        ]
        warp_mesh_by_planes_sdf!(mesh, plane_definitions, 0.3)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_6-plane_cut.vtu"))
        @test check_watertight(mesh).open_edges == 0
        @test count_inverted_exact(mesh) == 0
    end
end
