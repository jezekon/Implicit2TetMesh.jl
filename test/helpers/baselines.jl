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
const BEAM_NODES = 22385          # length(mesh.X)
const BEAM_TETS = 92675           # length(mesh.IEN)

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
const BEAM_BISECT_NODES = 22746   # length(mesh.X)
const BEAM_BISECT_TETS = 94640    # length(mesh.IEN)

const BEAM_DIH_BISECT = (
    min = 14.607,                 # smallest dihedral angle
    max = 148.575,                # largest dihedral angle
    lt5 = 0,                      # dihedral angles < 5 deg
    lt10 = 0,                     # dihedral angles < 10 deg
    gt140 = 61,                   # dihedral angles > 140 deg
    inverted = 0,                 # tets with non-positive float signed volume
)

# ------------------------------------------------------------------------------
# Beam, no planes, with the Etapa-10 relaxation post-pass (cut_points = :linear).
# Topology is unchanged (a relaxation only MOVES nodes), so the node/tet counts equal
# the no-relax beam; only the dihedral histogram moves. The gate guarantees min dihedral
# does not drop AND max dihedral does not grow vs the no-relax BEAM_DIH (min 11.276, max
# 153.211, gt140 81): both modes improve both tails. Re-freezable like the others.
# ------------------------------------------------------------------------------
# mode = :uniform (DistMesh size-equalizing springs): the surface re-projects onto the
# trilinear zero and sizes equalize, clearing every cap (gt140 81 -> 0).
const BEAM_RELAX_UNIFORM_DIH = (
    min = 18.373,                 # smallest dihedral angle (was 11.276)
    max = 131.616,                # largest dihedral angle  (was 153.211)
    lt5 = 0,
    lt10 = 0,
    gt140 = 0,                    # caps eliminated (was 81)
    inverted = 0,
)

# mode = :quality (quartet-style maximin smoothing): conservative two-sided gate lifts
# the worst angles a little and shrinks the cap tail (gt140 81 -> 9).
const BEAM_RELAX_QUALITY_DIH = (
    min = 12.883,                 # smallest dihedral angle (was 11.276)
    max = 148.378,                # largest dihedral angle  (was 153.211)
    lt5 = 0,
    lt10 = 0,
    gt140 = 9,                    # cap tail shrunk (was 81)
    inverted = 0,
)

# mode = :uniform, sizing = :curvature (mild element concentration where the surface
# bends; Persson 2006 gradient-limited h-field). Close to plain :uniform on the beam.
const BEAM_RELAX_CURVATURE_DIH = (
    min = 18.468,
    max = 132.999,
    lt5 = 0,
    lt10 = 0,
    gt140 = 0,
    inverted = 0,
)

# ------------------------------------------------------------------------------
# Beam, two cutting planes (Square(30) @ x=0, Square(5) @ x=60, warp_param 0.3).
# The plane warp only moves nodes, so the counts equal the no-planes mesh.
# ------------------------------------------------------------------------------
const BEAM_CUT_NODES = 22385      # length(mesh.X) after the plane cut
const BEAM_CUT_TETS = 92675       # length(mesh.IEN) after the plane cut

# ------------------------------------------------------------------------------
# Gripper, no cutting planes (opt-in: set I2TM_TEST_GRIPPER=1)
# ------------------------------------------------------------------------------
const GRIPPER_NODES = 414072      # length(mesh.X)
const GRIPPER_TETS = 1984321      # length(mesh.IEN)

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

# ------------------------------------------------------------------------------
# Etapa 11 -- convex boundary cap recovery (recover_boundary_caps!, default OFF).
# A cap is a 1:3 split of a boundary tet that under-cuts a convex bulge: it inserts ONE
# new node on the true zero level set and turns ONE tet into THREE, so the deltas versus
# the no-caps mesh always satisfy (tets - tets_off) == 2 * (nodes - nodes_off) and the
# boundary-face count rises by 2 per cap. Caps recover convex VOLUME (total volume rises
# toward the field's) but at a deliberate QUALITY cost -- the cap tets are slivers, which
# is why the pass is OFF by default and `:quality`/`:uniform` relaxation is the documented
# mitigation (see README "Convex cap recovery"). The pass is deterministic, so these counts
# are a sharp regression signal; re-freezable like the others.
# ------------------------------------------------------------------------------
# Beam, default CapRecoveryOptions (sagitta_frac 0.2). 643 caps: +643 nodes, +1286 tets.
const BEAM_CAPS_NODES = 23028     # length(mesh.X) with recover_caps on (was BEAM_NODES 22385)
const BEAM_CAPS_TETS = 93961      # length(mesh.IEN) with recover_caps on (was BEAM_TETS 92675)

# Unstructured sphere with sagitta_frac 0.05 (the default 0.2 barely fires on this finely
# resolved :bisection sphere -- only ~5 caps; 0.05 exercises the pass with 324 caps). Caps
# never worsen the geometric fidelity to the analytic sphere (new apices sit on the same
# interpolated zero set as the warped boundary), so SPHERE_U_GEO_FIDELITY still bounds it.
const SPHERE_U_CAPS_NODES = 20018  # length(mesh.X), caps on (was SPHERE_U_NODES 19694)
const SPHERE_U_CAPS_TETS = 100186  # length(mesh.IEN), caps on (was SPHERE_U_TETS 99538)
