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
