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
        @test opts.cut_points == :linear

        # Rejects an unsupported scheme, a negative warp_param, and an unknown cut-point mode.
        @test_throws ErrorException MeshGenerationOptions(scheme = "BCC")
        @test_throws ErrorException MeshGenerationOptions(warp_param = -0.1)
        @test_throws ErrorException MeshGenerationOptions(cut_points = :newton)

        # A valid non-default configuration is accepted.
        ok = MeshGenerationOptions(warp_param = 0.5, quality_export = true)
        @test ok.warp_param == 0.5
        @test ok.quality_export == true
        @test MeshGenerationOptions(cut_points = :bisection).cut_points == :bisection

        # Relaxation is OFF by default; a RelaxOptions wires in cleanly.
        @test MeshGenerationOptions().relax === nothing
        @test MeshGenerationOptions(relax = RelaxOptions()).relax isa RelaxOptions
    end

    @testset "RelaxOptions validation" begin
        # Defaults are as documented.
        r = RelaxOptions()
        @test r.mode == :uniform
        @test r.sizing == :uniform
        @test r.max_sweeps == 10
        @test r.band == 2
        @test r.omega == 0.3
        @test r.frozen == Int[]

        # Rejects unknown modes / sizings and non-positive numeric parameters.
        @test_throws ErrorException RelaxOptions(mode = :both)
        @test_throws ErrorException RelaxOptions(sizing = :foo)
        @test_throws ErrorException RelaxOptions(max_sweeps = 0)
        @test_throws ErrorException RelaxOptions(band = -1)
        @test_throws ErrorException RelaxOptions(omega = 0.0)
        @test_throws ErrorException RelaxOptions(alpha = -1.0)

        # A valid quality configuration with a freeze-set is accepted.
        q = RelaxOptions(mode = :quality, band = 1, frozen = [1, 2, 3])
        @test q.mode == :quality
        @test q.band == 1
        @test q.frozen == [1, 2, 3]
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

    # The :bisection cut-point mode (Labelle & Shewchuk §3.1) exists for input fields that
    # are NOT distance-like -- the beam's RBF-smoothed SDF is exactly such a field, so it
    # doubles as the regression input. The numbers differ from the :linear default BY
    # DESIGN (frozen separately in BEAM_DIH_BISECT); the mode's defining guarantee is
    # surface fidelity: every boundary vertex must sit on the trilinear zero set. With
    # :linear the beam's boundary strays up to ~0.13 field units from it (the
    # dented-flat-walls artifact).
    @testset "Beam pipeline (bisection cut points)" begin
        outdir = ensure_output_dir()
        mesh = generate_tetrahedral_mesh(
            beam_grid_file(),
            beam_sdf_file(),
            joinpath(outdir, "beam_bisection");
            options = MeshGenerationOptions(cut_points = :bisection),
        )

        # Invariants
        @test check_watertight(mesh).open_edges == 0
        @test count_inverted_exact(mesh) == 0
        @test count_components(mesh) == 1
        @test node_sdf_consistency(mesh).pristine_within_tol
        @test min_signed_volume(mesh) > 0

        # Surface fidelity: the boundary tracks the implicit surface (the reason this mode
        # exists). 1e-6 is generous -- 40 bisection halvings land within ~1e-12 of the zero.
        @test boundary_max_abs_sdf(mesh) <= 1e-6

        # Baselines (shared with the per-stage diagnostic's Stage 5)
        @test length(mesh.X) == BEAM_BISECT_NODES
        @test length(mesh.IEN) == BEAM_BISECT_TETS
        check_dihedral_baseline(dihedral_stats(mesh), BEAM_DIH_BISECT)
    end

    # Per-stage beam diagnostic (defines its own "Beam per-stage diagnostics"
    # testset and exports a VTU after each phase to test/output/).
    include(joinpath(@__DIR__, "test_beam_stages.jl"))

    # Etapa 8 -- pluggable SDF sources: unit checks, the beam round-trip through the
    # unstructured HEX8 path, and a genuinely-unstructured analytic sphere.
    include(joinpath(@__DIR__, "test_sdf_sources.jl"))

    # Etapa 10 -- optional gated relaxation. The gate is a CONSTRUCTION: an accepted move
    # never worsens either dihedral tail, so "min dihedral does not drop" and "max dihedral
    # does not grow" are hard, testable guarantees here. Build the pre-relax beam once and
    # relax deepcopies (the pre-relax mesh IS the frozen no-relax BEAM_DIH baseline).
    @testset "Relaxation (Etapa 10)" begin
        outdir = ensure_output_dir()
        base = generate_tetrahedral_mesh(beam_grid_file(), beam_sdf_file(),
                                         joinpath(outdir, "relax_base"))
        base_cv = surface_edge_cv(base)

        @testset "mode = :uniform" begin
            m = deepcopy(base)
            relax_mesh!(m, RelaxOptions(mode = :uniform))

            # Topology is unchanged by a relaxation.
            @test length(m.X) == BEAM_NODES
            @test length(m.IEN) == BEAM_TETS

            # Invariants.
            @test check_watertight(m).open_edges == 0
            @test count_inverted_exact(m) == 0
            @test count_components(m) == 1
            @test max_interior_sdf(m) < 0.0          # every interior vertex strictly inside

            # Hard gate guarantees vs the no-relax baseline (both tails improve).
            s = dihedral_stats(m)
            @test s.min >= BEAM_DIH.min               # min dihedral does not drop
            @test s.max <= BEAM_DIH.max               # max dihedral does not grow
            @test s.gt140 <= BEAM_DIH.gt140           # cap count does not grow
            @test s.inverted == 0

            # The point of :uniform -- element sizes equalize near the boundary.
            @test surface_edge_cv(m) < base_cv

            # Frozen histogram (re-freezable).
            check_dihedral_baseline(s, BEAM_RELAX_UNIFORM_DIH)

            # The pipeline `relax` option must give exactly the direct relax_mesh! result.
            mw = generate_tetrahedral_mesh(beam_grid_file(), beam_sdf_file(),
                joinpath(outdir, "relax_wired");
                options = MeshGenerationOptions(relax = RelaxOptions(mode = :uniform)))
            @test mw.X == m.X
            @test mw.IEN == m.IEN
        end

        @testset "mode = :quality" begin
            m = deepcopy(base)
            relax_mesh!(m, RelaxOptions(mode = :quality))

            @test length(m.X) == BEAM_NODES
            @test length(m.IEN) == BEAM_TETS
            @test check_watertight(m).open_edges == 0
            @test count_inverted_exact(m) == 0
            @test count_components(m) == 1
            @test max_interior_sdf(m) < 0.0

            s = dihedral_stats(m)
            @test s.min >= BEAM_DIH.min               # min dihedral strictly lifted here
            @test s.max <= BEAM_DIH.max
            @test s.gt140 <= BEAM_DIH.gt140           # >140 count does not grow
            @test s.inverted == 0
            check_dihedral_baseline(s, BEAM_RELAX_QUALITY_DIH)
        end

        @testset "sizing = :curvature" begin
            m = deepcopy(base)
            relax_mesh!(m, RelaxOptions(mode = :uniform, sizing = :curvature))

            @test check_watertight(m).open_edges == 0
            @test count_inverted_exact(m) == 0
            @test count_components(m) == 1
            @test max_interior_sdf(m) < 0.0

            s = dihedral_stats(m)
            @test s.min >= BEAM_DIH.min
            @test s.max <= BEAM_DIH.max
            @test s.inverted == 0
            check_dihedral_baseline(s, BEAM_RELAX_CURVATURE_DIH)
        end

        # Surface vertices stay EXACTLY on phi = 0. This is clean only when the pre-relax
        # surface is already on the trilinear zero, i.e. with cut_points = :bisection (the
        # default :linear puts surface nodes on linear cut points, off the zero by design,
        # and the gate leaves the un-moved ones there). Re-projection keeps them on it.
        @testset "surface stays on phi = 0 (:bisection)" begin
            bm = generate_tetrahedral_mesh(beam_grid_file(), beam_sdf_file(),
                joinpath(outdir, "relax_bisect");
                options = MeshGenerationOptions(cut_points = :bisection,
                                                relax = RelaxOptions(mode = :uniform)))
            @test check_watertight(bm).open_edges == 0
            @test count_inverted_exact(bm) == 0
            @test boundary_max_abs_sdf(bm) <= 1e-6    # every boundary vertex on the zero set
        end

        # Determinism: same input -> identical output (so relaxed baselines can be frozen).
        @testset "determinism" begin
            for mode in (:uniform, :quality)
                d1 = deepcopy(base); relax_mesh!(d1, RelaxOptions(mode = mode))
                d2 = deepcopy(base); relax_mesh!(d2, RelaxOptions(mode = mode))
                @test d1.X == d2.X
                @test d1.IEN == d2.IEN
            end
        end

        # Plane alignment runs AFTER relaxation, so cutting planes still come out applied.
        @testset "with cutting planes" begin
            planes = [
                PlaneDefinition([-1.0, 0.0, 0.0], [0.0, 10.0, 0.0], Square(30.0)),
                PlaneDefinition([1.0, 0.0, 0.0], [60.0, 2.0, 2.0], Square(5.0)),
            ]
            pm = generate_tetrahedral_mesh(beam_grid_file(), beam_sdf_file(),
                joinpath(outdir, "relax_planes");
                options = MeshGenerationOptions(warp_param = 0.3, plane_definitions = planes,
                                                relax = RelaxOptions(mode = :uniform)))
            @test check_watertight(pm).open_edges == 0
            @test count_inverted_exact(pm) == 0
            @test length(pm.IEN) == BEAM_CUT_TETS     # topology preserved through relax + cut
        end
    end

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

            # Etapa 10 at scale: the gripper's bad tails are the relaxation target. The
            # exact histogram is not frozen here (sensitive at 2M tets); the hard gate
            # guarantees + invariants are the meaningful checks. :uniform clears the
            # sub-10 deg slivers and shrinks the cap tail.
            rm = deepcopy(mesh)
            relax_mesh!(rm, RelaxOptions(mode = :uniform))
            s0 = dihedral_stats(mesh)
            s = dihedral_stats(rm)
            @test length(rm.X) == GRIPPER_NODES       # topology unchanged
            @test length(rm.IEN) == GRIPPER_TETS
            @test check_watertight(rm).open_edges == 0
            @test count_inverted_exact(rm) == 0
            @test count_components(rm) == 1
            @test max_interior_sdf(rm) < 0.0
            @test s.min >= s0.min                     # min dihedral does not drop
            @test s.max <= s0.max                     # max dihedral does not grow
            @test s.gt140 <= s0.gt140                 # cap count does not grow
            @test s.lt10 <= s0.lt10                   # sliver count does not grow
            @test surface_edge_cv(rm) < surface_edge_cv(mesh)
        end
    end
end
