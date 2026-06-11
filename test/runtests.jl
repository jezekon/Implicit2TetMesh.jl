# ==============================================================================
# Implicit2TetMesh.jl test suite
# ==============================================================================
# Two kinds of tests:
#   (a) INVARIANT  -- properties every correct output must satisfy (watertight
#       boundary, zero inverted tets, single component, node-SDF consistency,
#       positive volumes). These survive future pipeline changes unchanged.
#   (b) BASELINE   -- exact frozen numbers (node/tet counts, dihedral histogram)
#       in test/helpers/baselines.jl. The pipeline is deterministic, so exact
#       equality is safe; an INTENTIONAL change re-freezes baselines.jl.
#
# Run:                julia --project=. -e 'using Pkg; Pkg.test()'
# With gripper:       I2TM_TEST_GRIPPER=1 julia --project=. -e 'using Pkg; Pkg.test()'
# VTU artifacts land in test/output/ (git-ignored).

using Test
using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh
using Implicit2TetMesh.Modification
using Implicit2TetMesh.Utils

include(joinpath(@__DIR__, "helpers", "helpers.jl"))
include(joinpath(@__DIR__, "helpers", "baselines.jl"))

# Assert a dihedral histogram against a frozen baseline NamedTuple. Min/max are
# compared rounded to 3 decimals; the counts exactly.
function check_dihedral_baseline(stats, base)
    @test round(stats.min; digits = 3) == base.min
    @test round(stats.max; digits = 3) == base.max
    @test stats.lt5 == base.lt5
    @test stats.lt10 == base.lt10
    @test stats.gt140 == base.gt140
    @test stats.inverted == base.inverted
end

@testset "Implicit2TetMesh.jl" begin

    @testset "Options validation" begin
        # Defaults are as documented.
        opts = MeshGenerationOptions()
        @test opts.scheme == "A15"
        @test opts.warp_param == 0.3
        @test opts.plane_definitions === nothing
        @test opts.quality_export == false

        # Rejects an unsupported scheme and a negative warp_param.
        @test_throws ErrorException MeshGenerationOptions(scheme = "BCC")
        @test_throws ErrorException MeshGenerationOptions(warp_param = -0.1)

        # A valid non-default configuration is accepted.
        ok = MeshGenerationOptions(warp_param = 0.5, quality_export = true)
        @test ok.warp_param == 0.5
        @test ok.quality_export == true
    end

    @testset "Beam pipeline (no planes)" begin
        outdir = ensure_output_dir()
        mesh = generate_tetrahedral_mesh(
            beam_grid_file(),
            beam_sdf_file(),
            joinpath(outdir, "beam_no_planes"),
        )

        # Invariants
        @test check_watertight(mesh).open_edges == 0
        @test count_inverted_exact(mesh) == 0
        @test count_components(mesh) == 1
        @test node_sdf_consistency(mesh).pristine_within_tol
        @test min_signed_volume(mesh) > 0

        # Baselines
        @test length(mesh.X) == BEAM_NODES
        @test length(mesh.IEN) == BEAM_TETS
        check_dihedral_baseline(dihedral_stats(mesh), BEAM_DIH)
    end

    @testset "Beam pipeline (with cutting planes)" begin
        outdir = ensure_output_dir()
        plane_definitions = [
            PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
            PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0)),
        ]
        options = MeshGenerationOptions(
            warp_param = 0.3,
            plane_definitions = plane_definitions,
        )
        mesh = generate_tetrahedral_mesh(
            beam_grid_file(),
            beam_sdf_file(),
            joinpath(outdir, "beam_planes");
            options = options,
        )

        # Invariants after the cut
        @test check_watertight(mesh).open_edges == 0
        @test count_inverted_exact(mesh) == 0

        # Baselines for the cut mesh
        @test length(mesh.X) == BEAM_CUT_NODES
        @test length(mesh.IEN) == BEAM_CUT_TETS
    end

    # Per-stage beam diagnostic (defines its own "Beam per-stage diagnostics"
    # testset and exports a VTU after each phase to test/output/).
    include(joinpath(@__DIR__, "test_beam_stages.jl"))

    @testset "Gripper (opt-in)" begin
        if get(ENV, "I2TM_TEST_GRIPPER", "0") != "1"
            @info "Gripper suite skipped -- set I2TM_TEST_GRIPPER=1 to run it."
        else
            outdir = ensure_output_dir()
            mesh = generate_tetrahedral_mesh(
                gripper_grid_file(),
                gripper_sdf_file(),
                joinpath(outdir, "gripper_no_planes"),
            )

            # Invariants
            @test check_watertight(mesh).open_edges == 0
            @test count_inverted_exact(mesh) == 0
            @test count_components(mesh) == 1
            @test node_sdf_consistency(mesh).pristine_within_tol
            @test min_signed_volume(mesh) > 0

            # Baselines
            @test length(mesh.X) == GRIPPER_NODES
            @test length(mesh.IEN) == GRIPPER_TETS
            check_dihedral_baseline(dihedral_stats(mesh), GRIPPER_DIH)
        end
    end
end
