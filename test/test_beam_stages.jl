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
# The stages run with cut_points = :bisection (owner's choice): the beam's input
# is an RBF-smoothed, NON-distance field, and the per-stage VTUs are meant for
# visual inspection -- with :linear the flat walls come out dented (see
# Stencils.jl cut_edge!), which is exactly the artifact this diagnostic would be
# used to chase. The optional cap-recovery pre-pass (recover_caps!, default-off in
# the pipeline) is ALSO enabled here, so the inspected beam includes the recovered
# caps. Stage 5 therefore asserts the BEAM_DIH_BISECT_CAPS baseline, not the cap-free
# :bisection BEAM_DIH_BISECT used by the core no-planes pipeline test.
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
        warp!(mesh, "A15"; cut_points = :bisection)
        update_connectivity!(mesh; build_ine = false)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_2-warped.vtu"))
        @test all_finite_coords(mesh)
        nsc = node_sdf_consistency(mesh)
        @test nsc.pristine_within_tol     # untouched lattice nodes still match eval_sdf exactly
        @test nsc.n_on_surface > 0        # warp snapped some nodes onto the surface (node_sdf = 0)
    end

    # --- Stage 3: trim spikes / slice crossing tets --------------------------
    @testset "Stage 3 slice_ambiguous_tetrahedra!" begin
        # Cap recovery runs between warp's update_connectivity! and the slice, exactly where
        # the pipeline calls it (default-off there; always on in this diagnostic). It snaps
        # inward "+000 spike" apexes onto phi = 0 so the slice keeps them as quad-zero caps.
        recover_caps!(mesh)
        slice_ambiguous_tetrahedra!(mesh, "A15"; cut_points = :bisection)
        update_connectivity!(mesh; build_ine = false)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_3-sliced.vtu"))
        @test check_watertight(mesh).open_edges == 0
        @test all_indices_in_bounds(mesh)
        # Surface fidelity of the bisection mode: every boundary vertex sits on the
        # trilinear zero set (with :linear the beam strays up to ~0.13 -> dented walls).
        @test boundary_max_abs_sdf(mesh) <= 1e-6
    end

    # --- Stage 4: drop inverted elements -------------------------------------
    @testset "Stage 4 remove_inverted_elements!" begin
        remove_inverted_elements!(mesh)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_4-inverted_removed.vtu"))
        @test count_inverted_exact(mesh) == 0
    end

    # --- Stage 5: keep the largest component (final refresh builds INE) ------
    # With cap recovery on, this no longer equals the cap-free no-planes pipeline, so it
    # asserts its OWN baseline BEAM_DIH_BISECT_CAPS (the shared cap-free BEAM_DIH_BISECT is
    # left for the core bisection test). NOTE the thin tail: snapping apexes distorts a few
    # neighbour tets below 10 deg (min 5.826, lt10 = 4) -- those are solid neighbours, not
    # caps, so the [10,140] cap filter does not remove them.
    @testset "Stage 5 remove_isolated_components!" begin
        remove_isolated_components!(mesh, keep_largest = true)
        update_connectivity!(mesh)
        export_mesh_vtu(mesh, joinpath(outdir, "Beam-stage_5-components.vtu"))
        @test count_components(mesh) == 1
        @test check_watertight(mesh).open_edges == 0
        s = dihedral_stats(mesh)
        @test round(s.min; digits = 3) == BEAM_DIH_BISECT_CAPS.min
        @test round(s.max; digits = 3) == BEAM_DIH_BISECT_CAPS.max
        @test s.lt5 == BEAM_DIH_BISECT_CAPS.lt5
        @test s.lt10 == BEAM_DIH_BISECT_CAPS.lt10
        @test s.gt140 == BEAM_DIH_BISECT_CAPS.gt140
        @test s.inverted == BEAM_DIH_BISECT_CAPS.inverted
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
