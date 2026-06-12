# ==============================================================================
# Frozen regression baselines for Implicit2TetMesh.jl
# ==============================================================================
# These are EXACT numbers captured from the current dev HEAD by running the
# pipeline once. The pipeline is deterministic (byte-identical reruns verified in
# Etapa 5), so exact equality is a safe, sharp regression check.
#
# Min/max dihedral angles are floats; the suite compares them ROUNDED to 3
# decimals (the precision recorded throughout the project history), so store the
# rounded value here. The integer counts are compared exactly.
#
# RE-FREEZING: after an INTENTIONAL pipeline change (e.g. Etapa 6) these numbers
# move on purpose. Re-freeze them by running the suite, copying the new values
# here, and committing THIS ONE FILE together with the pipeline change so the
# baseline update is a single, reviewable diff.

# ------------------------------------------------------------------------------
# Beam, no cutting planes (final pipeline output; also the per-stage Stage 5 result)
# ------------------------------------------------------------------------------
const BEAM_NODES = 22370          # length(mesh.X)
const BEAM_TETS = 92532           # length(mesh.IEN)

# Dihedral histogram (min/max in degrees rounded to 3 decimals; counts exact):
const BEAM_DIH = (
    min = 11.276,                 # smallest dihedral angle
    max = 153.211,                # largest dihedral angle
    lt5 = 0,                      # dihedral angles < 5 deg
    lt10 = 0,                     # dihedral angles < 10 deg
    gt140 = 81,                   # dihedral angles > 140 deg
    inverted = 0,                 # tets with non-positive float signed volume
)

# ------------------------------------------------------------------------------
# Beam, no cutting planes, cut_points = :bisection (Labelle & Shewchuk §3.1).
# The numbers differ from the :linear baseline BY DESIGN: on the beam's smoothed
# (non-distance) SDF the bisected cut points sit on the true trilinear zero set,
# so more crossing-edge cuts survive (more tets) and the sliver-producing
# misplaced cut vertices disappear (better dihedral extremes).
# ------------------------------------------------------------------------------
const BEAM_BISECT_NODES = 22742   # length(mesh.X)
const BEAM_BISECT_TETS = 94594    # length(mesh.IEN)

const BEAM_DIH_BISECT = (
    min = 14.607,                 # smallest dihedral angle
    max = 148.575,                # largest dihedral angle
    lt5 = 0,                      # dihedral angles < 5 deg
    lt10 = 0,                     # dihedral angles < 10 deg
    gt140 = 61,                   # dihedral angles > 140 deg
    inverted = 0,                 # tets with non-positive float signed volume
)

# ------------------------------------------------------------------------------
# Beam, two cutting planes (Square(30) @ x=0, Square(5) @ x=60, warp_param 0.3).
# The plane warp only moves nodes, so the counts equal the no-planes mesh.
# ------------------------------------------------------------------------------
const BEAM_CUT_NODES = 22370      # length(mesh.X) after the plane cut
const BEAM_CUT_TETS = 92532       # length(mesh.IEN) after the plane cut

# ------------------------------------------------------------------------------
# Gripper, no cutting planes (opt-in: set I2TM_TEST_GRIPPER=1)
# ------------------------------------------------------------------------------
const GRIPPER_NODES = 414072      # length(mesh.X)
const GRIPPER_TETS = 1984316      # length(mesh.IEN)

const GRIPPER_DIH = (
    min = 9.537,                  # smallest dihedral angle
    max = 155.359,                # largest dihedral angle
    lt5 = 0,                      # dihedral angles < 5 deg
    lt10 = 11,                    # dihedral angles < 10 deg
    gt140 = 1136,                 # dihedral angles > 140 deg
    inverted = 0,                 # tets with non-positive float signed volume
)

# ------------------------------------------------------------------------------
# Etapa 8 -- genuinely-unstructured analytic sphere. A graded + interior-jittered
# conforming HEX8 block on [0,4]^3 (seed 7) carries the analytic field phi = |x-c|-r
# with c=(2,2,2), r=1.2; it is meshed through the FE path on a structured lattice
# (dx=0.15, padding 2, cut_points=:bisection). There is NO quartet oracle for
# unstructured inputs, so these are frozen from a verified-good run and re-freezable
# like the others. Counts are a sharp regression signal; the quality FLOOR and the
# geometric fidelity to the TRUE sphere are the meaningful correctness checks.
# (|eval_sdf| on the boundary is NOT ~0 here -- unlike the structured paths -- because
# this graded+jittered hex block is not a perfect tiling: warped HEX8 faces overlap by
# ~1e-5, so a boundary cut point sits strictly inside >1 element with different values.
# A property of the distorted input, not a defect; see test/test_sdf_sources.jl.)
# ------------------------------------------------------------------------------
const SPHERE_U_NODES = 19694      # length(mesh.X)
const SPHERE_U_TETS = 99538       # length(mesh.IEN)
const SPHERE_U_GEO_FIDELITY = 0.05  # max boundary-vertex distance to |x-c|=r (meas. 0.0352)
const SPHERE_U_VOL_RELERR = 0.07    # |tet volume - 4/3 pi r^3| / (4/3 pi r^3) (meas. 0.0516)
