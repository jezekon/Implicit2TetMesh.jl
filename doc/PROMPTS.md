# Stage prompts — systematic alignment of Implicit2TetMesh.jl with the published isosurface‑stuffing algorithm

## How to use this file
1. Work happens on git branch **`dev`**.
2. For each stage: open a fresh Claude Code session at the repo root, paste the
   **Common context** block first, then the **stage prompt**.
3. Do the stages **in order**. Review between stages.
4. Every stage must **start** by gathering context and confirming the direction is
   consistent with the literature, and **end** by running regression, updating
   `README.md`, committing to `dev`, and stopping for review.

**Status (updated 2026-06-11):**
  • Etapas 1-5 — DONE, committed on `dev` (Etapa 4 = 25c6f70, Etapa 5 = 0bc8409; awaiting review).
  • Volume correction — REMOVED (standalone cleanup, c224f05); the pipeline no longer has a
    correct_mesh_volume! step, so any "compare before correct_mesh_volume!" note below is moot.
  • Etapa 9 — only the prompt text lives in this file (12e9bc9); no code yet.
  • NEXT coding stage: **Etapa 6** (optional adaptive octree sizing; start of Phase 2).
    Etapas 6-9 are optional/later. NOTE: the "Key files" line numbers below have drifted since
    Etapa 5 (Stencils.jl / RemoveIsolatedComponents.jl grew) — re-grep before relying on them.

Etapa 0 (validation experiment) is also complete. It proved that replacing the
current Newton‑to‑isosurface warp with quartet‑style edge warping collapses the slivers
(min dihedral 0.001°→8.19°, angles <5°: 95→0, bubbles 35→4). Its numbers are the
baselines below.

---

## 🔁 Common context (paste at the top of EVERY stage session)

````text
PROJECT CONTEXT — Implicit2TetMesh.jl

Implicit2TetMesh.jl is a Julia package that builds tetrahedral meshes from an SDF using an
isosurface-stuffing-style "trim spikes" approach on an A15 acute lattice.

CONVENTION: phi < 0 = INSIDE, phi > 0 = OUTSIDE, phi = 0 = ON the surface — matching the
reference implementations (quartet / Labelle / isostuffer). NOTE: the project ORIGINALLY used
the OPPOSITE convention (positive = inside); Etapa 1 flips it. Etapas 2+ assume phi < 0 = inside.

Pipeline (src/TetMeshGenerator.jl -> generate_tetrahedral_mesh):
  generate_mesh! -> warp! -> update_connectivity! -> slice_ambiguous_tetrahedra!
  -> update_connectivity! -> remove_inverted_elements! -> remove_isolated_components!(keep_largest)
  -> update_connectivity!
  (slice_ambiguous_tetrahedra! lost its 3rd arg in Etapa 3; the optional correct_mesh_volume! tail
   was removed with volume correction — there is no post-meshing node-mover anymore.)

Key files (line numbers as of 2026-06-10, post-Etapa-4 + volume removal; re-grep if they drift):
  src/GenerateMesh/TetGenerator.jl    — process_cell_A15! (:11), generate_mesh! (:194), warp! (:250),
                                        connectivity: merge_duplicate_nodes! (:120), cleanup_unused_nodes!
                                        (:153), create_INE! (:181), update_connectivity! (:313)
  src/GenerateMesh/Stencils.jl        — slice_ambiguous_tetrahedra! (:21), cut_edge! (:69),
                                        apply_stencil_trim_spikes! (:179, trim cases :311-396, all-on-surface
                                        Case 3 deferred :202), Labelle §3.4 surface-candidate decision incl.
                                        centroid test (:436-512), orientation check_tetrahedron_orientation
                                        (:124) / fix_tetrahedron_orientation! (:152), remove_inverted_elements!
                                        (:539).  [experimental NZZZ / process_nzzz_case! / Schläfli: removed in Etapa 3]
  src/GenerateMesh/NewCases-Experimental.jl — create_warping_params (:54), compute_dihedral_angle_range (:76)
  src/GenerateMesh/Schemes/A15Scheme.jl — A15 tile: tile_ref (:3), tetra_connectivity (:35)
  src/Fundamentals/BlockMesh.jl       — mutable struct (:4) (X, IEN, INE, SDF, node_sdf, node_hash, grid_step, grid_tol)
  src/Fundamentals/SDFOperations.jl   — get_cell_sdf_values (:4), eval_sdf (:24), compute_gradient (:89/:108)
  src/Modification/RemoveIsolatedComponents.jl — face->elements map pattern (:34-49)  [reuse for adjacency]
  src/Modification/ModifyResultingMesh.jl      — find_surface_nodes (:123), surface_faces count==1 (:146),
                                                 is_on_plane (:2)

Reference materials (read the relevant parts per stage):
  Literature/quartet-original/src/make_tet_mesh.cpp — THE A15 reference (Bridson & Doran):
        warp_vertices (152-193), trim_spikes (221-346), remove_exterior_tets (350-366)
  Literature/2007_Labelle_IsosurfaceStuffing.pdf — algorithm + guarantees: §3.2 warp (p.4), §3.3-3.4 (p.5)
  Literature/isostuffer/src/ — BCC graded-octree reference: OctTree.*, IsoStuffer.hpp (stencil_match)

OVERARCHING GOAL: make the mesher CONFORM to the published quartet/Labelle algorithm — robust and
fast — as a foundation to build on.
ROOT CAUSE already proven (Etapa 0): the current warp! diverges from the published algorithm (it
warps every nearby node to the true isosurface via Newton, instead of snapping a vertex to the
linear cut point on an incident edge). That divergence causes the slivers and the bubbles.

DECISIONS ALREADY MADE:
  • A15 only (Schläfli removed in Etapa 3).
  • phi < 0 = inside, matching the originals (Etapa 1).
  • Exact predicates via ExactPredicates.jl (Etapa 5).
  • Graded refinement is an OPTIONAL feature, OFF by default; the uniform mesher stays the default
    path (Etapas 6-7).
  • The SDF source will become pluggable (Etapa 8): structured grid (trilinear, current path) AND
    unstructured conforming HEX8 mesh (FE shape-function interpolation), with SIMP/level-set input
    adapters. IMPORTANT for Etapas 2-5: never add new direct mesh.SDF[i,j,k] reads — access the
    field ONLY via eval_sdf / get_cell_sdf_values / node_sdf, so the interface seam stays narrow.

BEFORE STARTING (do this first, EVERY stage):
  1. Read the code files and the reference sections listed in the stage.
  2. Confirm the planned change matches the published algorithm (quartet/Labelle). Write 3-5 lines
     stating WHY the approach is correct and literature-consistent.
  3. If anything conflicts with the literature, or the current code differs from what the stage
     assumes, STOP and report before changing any code.

CODE STYLE (the maintainer is at an intermediate Julia level — readability first):
  • Match the existing style of the file you edit. Prefer clear, explicit, well-structured code
    over advanced/idiomatic Julia (avoid heavy metaprogramming, deep type parametrization, clever
    broadcasting). Use explicit loops and descriptive names. Every line must be readable.
  • All code comments and docstrings in ENGLISH. Add a short docstring to each new function.

CLOSEOUT (do this at the end, EVERY stage):
  1. Run the regression harness on beam AND gripper; compare to baseline / previous stage.
  2. Update README.md so it matches the CURRENT implementation (options, schemes, behavior).
  3. Commit to branch `dev` with a clear message (create `dev` from the current branch if needed).
     Do not push unless asked.
  4. Present a short report: what changed, regression table, README changes. STOP for review.
  Keep changes scoped to the current stage; do not start later stages.

REGRESSION HARNESS (save as scratch/measure.jl; run: julia --project=. scratch/measure.jl beam).
Adapt the pipeline calls to the CURRENT code (e.g. after Etapa 3 the slice loses its 3rd arg).
The dihedral statistics are sign-agnostic (purely geometric), so they are valid in any convention.
```julia
using Implicit2TetMesh, Implicit2TetMesh.Fundamentals, Implicit2TetMesh.GenerateMesh, Implicit2TetMesh.Modification
using JLD2, LinearAlgebra
nm = isempty(ARGS) ? "beam" : ARGS[1]
f = Dict("beam"=>("Z_beam_HEX8_FineGrid_B-1.0_smooth-1.jld2","Z_beam_HEX8_FineSDF_B-1.0_smooth-1.jld2"),
         "gripper"=>("Z_robot_gripper_HEX8_FineGrid_B-2.0085_smooth-1.jld2","Z_robot_gripper_HEX8_FineSDF_B-2.0085_smooth-1.jld2"))
gf,sf=f[nm]; @load joinpath("data",nm,gf) fine_grid; @load joinpath("data",nm,sf) fine_sdf
function dih(mesh)
  mn=180.0;mx=0.0;l5=0;l10=0;g140=0;inv=0
  for t in mesh.IEN
    v=[mesh.X[t[k]] for k=1:4]; dot(v[2]-v[1],cross(v[3]-v[1],v[4]-v[1]))<=0 && (inv+=1)
    F=[(1,2,3,4),(1,2,4,3),(1,3,4,2),(2,3,4,1)];N=[];ok=true
    for (a,b,c,o) in F; cr=cross(v[b]-v[a],v[c]-v[a]);n=norm(cr);n<1e-14&&(ok=false;break)
      u=cr/n;ct=(v[a]+v[b]+v[c])/3;dot(u,v[o]-ct)>0&&(u=-u);push!(N,u);end
    ok||continue
    for i=1:3,j=(i+1):4;a=acos(-clamp(dot(N[i],N[j]),-1,1))*180/pi
      mn=min(mn,a);mx=max(mx,a);a<5&&(l5+=1);a<10&&(l10+=1);a>140&&(g140+=1);end
  end; (mn,mx,l5,l10,g140,inv)
end
m=BlockMesh(fine_sdf,fine_grid)
generate_mesh!(m,"A15"); warp!(m,"A15"); update_connectivity!(m)
slice_ambiguous_tetrahedra!(m,"A15",false); update_connectivity!(m)
remove_inverted_elements!(m); remove_isolated_components!(m,keep_largest=true); update_connectivity!(m)
s=dih(m); println("$nm: tets=$(length(m.IEN)) min=$(round(s[1];digits=3)) max=$(round(s[2];digits=3)) <5=$(s[3]) <10=$(s[4]) >140=$(s[5]) inv=$(s[6])")
```

BASELINE (beam, A15):
  • Current repo (before this work): tets 94,524 | min 0.001° | <5° 95 | <10° 210 | >140° 413 | isolated 35 | inverted 0
  • Etapa-0 edge-warp prototype:     tets 92,675 | min 8.19° | <5° 0  | <10° 23  | >140° 273 | isolated 4  | inverted 0
Etapa 1 must reproduce the "current repo" numbers EXACTLY (output-preserving). Etapa 2 should reach
≈the edge-warp prototype. Etapas 3-4 should improve it further. Always also run gripper.

CURRENT HEAD of dev (after Etapa 4 + volume removal — the numbers Etapa 5 must preserve):
  • beam:    tets 92,532    | min 11.276° | <10° 0  | >140° 81   | inverted 0
  • gripper: tets 1,984,316 | min 9.537°  | <10° 11 | >140° 1136 | inverted 0
````

---

## 🧪 Optional cross-validation against the original quartet implementation

Not a stage of its own — run it as an extra check at the end of **Etapa 2** (warp port) and
again after **Etapa 4** (full trim pipeline). It is a MANUAL verification workflow, not CI:
it needs a C++ build and a slow run. Referenced from the VERIFY steps of Etapas 2-4 and 8.

````text
WHY THIS WORKS: quartet's entry point make_tet_mesh(mesh, sdf, optimize, ...) (make_tet_mesh.h)
takes the SDF grid directly (Array3f phi + origin + dx) — the same input our BlockMesh consumes.
So we can bypass quartet's main.cpp (which computes the SDF from an .obj) and feed quartet
EXACTLY the same SDF grid as our mesher, then compare two meshes built from identical input.

SCOPE / LIMITS (be explicit about these in any report):
  • Works ONLY for STRUCTURED-grid inputs (beam, gripper, synthetic SDFs). quartet cannot consume
    an unstructured HEX8 field, so for Etapa 8's unstructured inputs there is NO quartet oracle —
    the only check there is the quality histogram (see Etapa 8).
  • Run quartet with optimize=false (we have no optimization pass) and without feature matching.
  • Compare our FINAL mesh: volume correction was REMOVED, so there is no correct_mesh_volume!
    step to exclude any more (quartet likewise has none).
  • Do NOT expect vertex-by-vertex / connectivity identity: quartet sizes its grid and enumerates
    A15 tiles its own way. The right altitude is aggregate metrics + boundary-surface distance.
  • isostuffer is BCC + graded octree — a different lattice; it is NOT comparable to the A15
    output. It stays reference material for Etapas 6-7 only.

SETUP (one-time; keep everything under scratch/):
  1. scratch/export_sdf.jl — dump the IN-MEMORY SDF (the negated field, phi<0=inside, as stored
     by BlockMesh) plus origin and dx for beam/gripper into a simple binary file. Do NOT dump the
     raw on-disk field — it still has the old positive=inside sign.
  2. scratch/quartet_driver.cpp (~50 lines) — read the dump, build an SDF (sdf.h), call
     make_tet_mesh(mesh, sdf, /*optimize=*/false), write the .tet output. Compile against
     Literature/quartet-original/src using its existing Makefile or a one-line g++ command.
  3. scratch/compare_quartet.jl — load both meshes and print side-by-side:
       - the SAME dih() metrics as the regression harness (tet count, min/max dihedral,
         <5°, <10°, >140°, inverted) for both meshes;
       - total mesh volume of each vs the SDF volume;
       - boundary-surface distance: mean + Hausdorff distance from our boundary vertices to
         quartet's boundary triangles AND vice versa (both meshes approximate the same zero
         isosurface; expect agreement within ~dx/10).

WHAT TO CHECK, BY STAGE:
  • After Etapa 2 (the strongest test — the warp is a 1:1 port): temporarily instrument quartet
    (printf) to dump vertex positions after warp_vertices; on the same lattice our warp must
    produce IDENTICAL displacements up to float tolerance.
  • After Etapa 3: per-case trim counts (+++-/++--/+0--/...) should match quartet's counts on the
    same input (instrument trim_spikes with counters).
  • After Etapa 4: aggregate metrics close to quartet's, inverted = 0 in both, boundary surfaces
    coincident within tolerance.
  • Any large deviation = a porting bug; find it before moving to the next stage.

FUTURE AUTOMATED TESTS: once a cross-validated run produces trusted reference numbers for
beam/gripper (tet count, dihedral histogram, volume), freeze them as expected values in test/ so
regressions are caught automatically WITHOUT the C++ build. For unstructured-input cases (Etapa 8)
no quartet oracle exists — instead, once the first verified-good unstructured run is approved by
the owner, freeze ITS histogram/quality metrics as the automated regression baseline.
````

---

## Etapa 1 — Adopt the standard sign convention: phi < 0 = inside  ✅ DONE

````text
ETAPA 1 — Flip the project to phi < 0 = INSIDE (output-preserving refactor)

GOAL: Change the project's SDF convention from "positive = inside" to "phi < 0 = inside" so it
matches the published algorithm and the reference code (quartet/Labelle/isostuffer). This makes
all later ports direct (no sign flipping). It is OUTPUT-PRESERVING: the generated mesh must be
identical before and after (a sign flip is just a relabeling).

BEFORE STARTING: follow the "BEFORE STARTING" checklist in the common context. Confirm from
quartet (make_tet_mesh.cpp: vphi<=0 means inside/kept) and Labelle §3 that phi<0=inside is the
reference convention.

IMPLEMENT:
  • Negate the SDF once at the source so the in-memory field uses phi<0=inside. Cleanest place:
    the BlockMesh constructor (src/Fundamentals/BlockMesh.jl) — store mesh.SDF = -Float64.(fine_sdf).
    Do NOT modify the data files on disk.
  • Flip EVERY sign-dependent decision to the new convention. Audit at minimum:
      - TetGenerator.jl: process_cell_A15! keep-tests (x >= -tol, x >= 0); warp! pass selection
        (sdf>0 / sdf<0). IMPORTANT: warp! does "inside first, then outside" in two passes — swap
        the comparisons so it still processes inside (now sdf<0) first, to preserve behavior/order.
      - Stencils.jl: the is_*_neg / is_*_pos / is_*_zero classification, cut_edge! / interpolate_zero
        polarity checks, every case branch, and the inline ZZZZ centroid test (eval_sdf(centroid) > 0).
      - SDFOperations.jl: compute_gradient / Newton step direction if used to move "toward" the surface
        (note: a global sign flip leaves the Newton step f/|grad|^2 * grad invariant — verify).
      - Modification/: CalcVolumeFromSDF, CorrectMeshVolume (correct_mesh_volume! receives the RAW
        fine_sdf at TetMeshGenerator.jl:112 — make it consistent with the negated convention),
        boundary-plane logic, anything that uses the SDF sign.
      - ExportMesh.jl / quality export if it uses the sign.
  • Update all English comments/docstrings that state "positive = inside" etc.

VERIFY: the regression harness must produce IDENTICAL numbers to the "current repo" baseline on
beam (94,524 tets, min 0.001°, <10° 210, >140° 413, isolated 35, inverted 0) and unchanged gripper
output. Any difference means a sign site was missed — find it. Then do the CLOSEOUT (incl. README:
update the convention description).
````

---

## Etapa 2 — Warp realignment (edge-based warp to cut points)  ✅ DONE

````text
ETAPA 2 — Replace warp! with an edge-based warp to cut points (port of quartet warp_vertices)
(depends on Etapa 1; the project is now phi<0=inside, so this ports quartet directly with no sign flip)

GOAL: Replace the node-based Newton-to-isosurface warp with an EDGE-BASED warp that snaps a lattice
vertex onto the LINEAR cut point of an incident sign-crossing edge when that cut point is within a
threshold of the vertex. Faithful port of quartet warp_vertices (make_tet_mesh.cpp:152-193).
Etapa 0 proved this collapses slivers; land it properly now.

BEFORE STARTING: follow the common-context checklist. Read quartet warp_vertices (152-193), the
current warp! / warp_node_to_isocontour! (TetGenerator.jl:269-342), and Labelle §3.2 (p.4).

IMPLEMENT (direct port; phi<0=inside):
  • Iterate edges from mesh.IEN. For an edge (i,j) with opposite signs
    ((phi_i < 0 && phi_j > 0) || (phi_i > 0 && phi_j < 0)):
        alpha = phi_i / (phi_i - phi_j)        # fraction from i where phi = 0, in (0,1)
        if alpha < threshold:        warp i -> cut point, disp_i = alpha * (X_j - X_i)
        elseif alpha > 1 - threshold: warp j -> cut point, disp_j = (1 - alpha) * (X_i - X_j)
    Each vertex keeps its CLOSEST qualifying cut point (smallest alpha * |edge|^2).
    Apply all warps at once: X[v] += disp_v; node_sdf[v] = 0.0.
  • Move to the LINEAR cut point, NOT the true isosurface (no Newton).
  • threshold = 0.3 (quartet's value); expose it as a function parameter so Etapa 3 can tune it.
  • Keep the public signature warp!(mesh, scheme) so the pipeline/call site is unchanged.
  • Remove the now-unused Newton warp logic (warp_node_to_isocontour!) if nothing else needs it.
  • Write the code in clear explicit loops with English comments (see CODE STYLE).

VERIFY: harness on beam AND gripper. Beam should land near the Etapa-0 prototype (≈92,675 tets,
min ≈8.2°, <5° 0, <10° ≈23, >140° ≈273, isolated ≈4, 0 inverted). Report a before/after table.
RECOMMENDED EXTRA: run the quartet cross-validation (see the "🧪 Optional cross-validation"
section) — on the same lattice and SDF, the edge-warp displacements must match quartet's
warp_vertices exactly (up to float tolerance). Structured-grid inputs only.
CLOSEOUT (README: describe the new warping step).
````

---

## Etapa 3 — Trim/stencil alignment with quartet + remove Schläfli  ✅ DONE

````text
ETAPA 3 — Align slicing 1:1 with quartet trim_spikes; remove the experimental and Schläfli paths
(depends on Etapa 2)

GOAL: Make the cut tetrahedra well-shaped by matching quartet trim_spikes exactly; remove the
non-conforming experimental branches; remove the Schläfli scheme; tune the warp threshold.

BEFORE STARTING: follow the common-context checklist. Read quartet trim_spikes (221-346) +
remove_exterior_tets (350-366); the current apply_stencil_trim_spikes! (Stencils.jl:389-646). The
case enumeration already mirrors quartet (NNNP/NPPP/NNPP/NZPP/NNZP/NZZP <-> +++-/+---/++--/+0--/
++0-/+00-) and the inline ZZZZ (:457-468) mirrors remove_exterior_tets. (Reminder: project is now
phi<0=inside, so the Julia cases line up with quartet's signs directly.)

TASKS:
  1. Audit each case branch against quartet: vertex ordering, the QUAD-SPLIT diagonal choice
     (quartet splits each quad consistently via its sorted pqrs so face-adjacent tets stay
     compatible — verify the Julia split matches; this is critical for a WATERTIGHT boundary), and
     orientation handling. Fix any divergence.
  2. Remove the experimental NNNZ/NNZZ/NZZZ branches (Stencils.jl:506-537) and process_nzzz_case!
     (:294-380). With Etapa 2's warp, these "no interior vertex" cases are correctly discarded as
     exterior (= quartet's vphi[s]==0 branch). Remove the experimental_nzzz flag and its plumbing
     (TetMeshGenerator.jl options + Stencils.jl signatures). KEEP compute_dihedral_angle_range
     (Etapa 4 uses it). Thin features are handled later by OPTIONAL refinement (Etapas 6-7) — not here.
  3. Remove Schläfli entirely: process_cell_Schlafli! (TetGenerator.jl:118-161),
     Schemes/SchlafliScheme.jl, the scheme branch in generate_mesh! and warp!, and Schläfli handling
     in create_warping_params. Make A15 the only path; update tests/docstrings.
  4. Tune the warp threshold (from Etapa 2) toward the published values to tighten angles; report
     the trade-off (smaller threshold = fewer warps = more cuts; larger = more warps).

VERIFY: harness on beam + gripper; residual <10° / >140° counts should drop vs Etapa 2; the boundary
must be watertight (check the extracted surface / getBoundary for cracks). Before/after table.
OPTIONAL: cross-check the per-case trim counts against quartet on the same SDF (see the
"🧪 Optional cross-validation" section).
CLOSEOUT (README: remove Schläfli and experimental_nzzz from the documented options).
````

---

## Etapa 4 — Quadruple-zero / bubble prevention (Labelle §3.4)  ✅ DONE (awaiting review)

````text
ETAPA 4 — Principled boundary handling per Labelle §3.4 (depends on Etapas 2-3)

GOAL: Upgrade the surface-tet ("quadruple-zero") handling from a centroid-only test to the full
published heuristic, eliminating bubbles the principled way (replacing the coarse global
remove_isolated_components! as the PRIMARY mechanism; keep it only as a final safety net).

BEFORE STARTING: follow the common-context checklist. Read Labelle §3.3-3.4 (p.5) — the "four
options" and the rule: discard a surface-only tet if it is inverted or its dihedral angles are
poor; of the nicely-shaped ones, RETAIN if all four faces adjoin output tets, DISCARD if none do
(prevents "bubbles"), decide the rest by the SDF sign at the centroid. Read quartet
remove_exterior_tets (350-366). Reuse the face->elements map (RemoveIsolatedComponents.jl:31-49)
and boundary-face extraction (ModifyResultingMesh.jl:122-151).

IMPLEMENT:
  • Build a face->count map of the interior mesh produced by the trim step. For each surface-only
    candidate tet, apply in order:
        inverted OR min-dihedral < bound OR max-dihedral > bound  -> discard
        all 4 faces adjoin interior tets                          -> retain
        0 faces adjoin                                            -> discard (bubble)
        else                                                      -> centroid SDF test (inside -> keep)
    Use compute_dihedral_angle_range for the dihedral test (bounds from create_warping_params).
  • If a candidate cannot see its neighbours in a single pass, use a deferred/two-pass structure:
    collect candidates, build the map from the solid mesh, then decide.
  • Keep remove_isolated_components! as a final safety net (it should now find ≈nothing).
  • Cleanups from the earlier audit: add any new node to the mesh only AFTER all checks pass (no
    orphan nodes); remove dead parameters.
  • Clear loops + English comments (CODE STYLE).

VERIFY: harness on beam + gripper. Target: isolated/bubble count -> 0 (beam was 4 after Etapa 2),
no quality regression, watertight. Before/after table.
RECOMMENDED EXTRA: run the full quartet cross-validation (see the "🧪 Optional cross-validation"
section) — same SDF into quartet (optimize=false), compare aggregate metrics, volume and
boundary-surface distance; compare our mesh BEFORE correct_mesh_volume!. If the numbers agree,
freeze them as automated regression values in test/.
CLOSEOUT (README: describe the boundary handling).
````

---

## Etapa 5 — Exact predicates + performance  ✅ DONE (awaiting review)

````text
ETAPA 5 — Robust geometric predicates + speed (depends on Etapas 2-4)

GOAL: Base orientation/inversion decisions on EXACT predicates (no tolerance guessing) and remove
redundant work.

CURRENT STATE (verified 2026-06-10): Etapas 1-4 + volume removal are on dev, tree clean.
ExactPredicates is ALREADY resolved in Manifest.toml (transitive), but NOT in Project.toml [deps] —
add it as a DIRECT dependency.

BEFORE STARTING: follow the common-context checklist. Read the current float tests (all in
src/GenerateMesh/Stencils.jl): check_tetrahedron_orientation (:124, returns
dot(a,cross(b,c)) > 1e-12 at :143), fix_tetrahedron_orientation! (:152), the inline orientation
check inside the slice (:488), and remove_inverted_elements! (:539, its own det computations at
:570/:594/:616/:637). quartet uses Shewchuk's exact predicates (predicates.cpp).

TASKS:
  1. Add ExactPredicates to Project.toml [deps] (it is already resolved in the Manifest). VERIFY
     its 3D orientation API (e.g. orient(a,b,c,d) returning -1/0/1), confirm the sign convention
     matches the current code (positive determinant = correct orientation), and replace the float
     SIGN tests for orientation/inversion with it. UNIFY: make ONE exact helper (e.g.
     is_positively_oriented(mesh, tet)) and call it from ALL the sites above, so the sign decision
     is not duplicated across five places. Keep a tiny tolerance ONLY for near-zero-volume
     degeneracy removal, but decide the SIGN exactly.
  2. Performance: the pipeline calls update_connectivity! 3x (TetGenerator.jl:313 =
     cleanup_unused_nodes!(:153) + merge_duplicate_nodes!(:120) + create_INE!(:181); pipeline call
     sites TetMeshGenerator.jl:81/:85/:90). Identify which rebuilds are actually required and
     collapse the rest. Reduce repeated eval_sdf where node_sdf is already cached. Time mesh
     generation on gripper (@time / a simple timer) before & after; report the speedup.
  3. Guarantee 0 inverted tets via the exact predicate (not tolerance).
  Keep wrappers small and readable (CODE STYLE); the maintainer should follow how the predicate is called.
  Common-context seam reminder (Etapa 8): do NOT add new direct mesh.SDF reads — access the field
  only via eval_sdf / get_cell_sdf_values / node_sdf.

VERIFY: scratch/measure.jl on beam + gripper — numbers must match the CURRENT HEAD baseline within
rounding, inverted = 0 EXACTLY:
    beam:    tets 92,532    | min 11.276° | <10° 0  | >140° 81   | inverted 0
    gripper: tets 1,984,316 | min 9.537°  | <10° 11 | >140° 1136 | inverted 0
Also run: cd test && julia --project=.. runtests.jl. Report timing before/after.
CLOSEOUT (README: note exact predicates + the new ExactPredicates dependency).
````

---

## Etapa 6 — Adaptive octree sizing (OPTIONAL feature; start of Phase 2)

````text
ETAPA 6 — Optional adaptive octree sizing for thin features (depends on a solid base, Etapas 1-5)

GOAL: Add an OPTIONAL mode that detects where the uniform A15 lattice under-resolves the geometry
(thin walls, high curvature) and produces a graded octree / target-level field that refines there.
This is the publication-backed way to capture thin features (Labelle §3.1 & §6: use a finer lattice
where the surface is not resolved). The uniform mesher remains the DEFAULT; refinement is opt-in.

BEFORE STARTING: follow the common-context checklist. Read Labelle §3.1 (p.3) and §6 (adaptivity);
isostuffer OctTree.h/.cpp (multi-level grid, level 0 = finest). Confirm the criterion you pick is
consistent with the literature's "feature size vs grid spacing" reasoning.

TASKS:
  • Add an option to MeshGenerationOptions to enable refinement (default OFF), e.g. a small
    RefinementOptions struct (max levels, criterion threshold). When OFF, behaviour is exactly the
    uniform pipeline.
  • Define an under-resolution criterion from the SDF (start simple/measurable): local feature size
    / distance-to-medial-axis (opposite-sign surface within ~1 cell, or |grad(SDF)| drop), or "a
    cell where the uniform pass would lose interior material / produce a quadruple-zero".
  • Build an octree (or per-cell target level) refining near such cells, 2:1 balanced (<=1 level
    difference across face-neighbours) so Etapa 7 transitions stay tractable.
  • Output a structure the lattice fill can consume (cell -> level). Do NOT change the stencils yet.
  • Add a synthetic thin-plate SDF (thinner than dx) as a demonstrable test; gripper is the real stress.

VERIFY: with refinement OFF, beam/gripper output is unchanged from Etapa 5. With it ON, thin regions
are refined and the tree is 2:1 balanced. CLOSEOUT (README: document the optional refinement and that
it is OFF by default).
````

---

## Etapa 7 — Graded fill + transition stencils (OPTIONAL; Phase 2)

````text
ETAPA 7 — Graded fill + level-transition stencils (OPTIONAL refinement; depends on Etapa 6)

⚠️ DECISION TO RESOLVE BEFORE CODING (present options to the owner):
Graded isosurface stuffing with transition stencils is FULLY PUBLISHED only for the BCC lattice
(Labelle §6; isostuffer implements it), NOT for A15 (quartet is uniform-only; Doran 2013 is a
1-page poster without A15 transition stencils). Two paths:
  (a) Design A15 transition stencils ourselves — keeps "A15 only", but research-grade.
  (b) Use BCC for the graded mode (proven, isostuffer-aligned) — accepts a second tile type, used
      only when refinement is enabled.
RECOMMENDED: first prototype grading on the proven BCC path (port isostuffer's level-transition
handling) to validate the machinery end-to-end, THEN decide whether A15 transitions are worth it.
Get the decision before committing.

BEFORE STARTING: follow the common-context checklist. Read isostuffer IsoStuffer.hpp (octree fill +
stencil_match at level boundaries) + tables.h; Labelle §6.

TASKS:
  • This is part of the OPTIONAL refinement mode (default OFF). Fill the finest cells with the chosen
    tile and generate transition elements between levels so the mesh stays CONFORMING (shared faces
    match exactly — no T-junctions / cracks).
  • Reuse the same warp (Etapa 2), boundary handling (Etapa 4), and exact predicates (Etapa 5) on the
    graded mesh.
  • Keep the code readable and English-commented (CODE STYLE); transition logic is the hardest part —
    add diagrams/comments explaining the level-boundary cases.

VERIFY: the synthetic thin-feature case is now captured (it was lost in the uniform mesher);
watertight, 0 inverted, dihedral within bounds; beam unchanged in uniform regions and with
refinement OFF. CLOSEOUT (README: document the graded mode and the chosen tile).
````

---

## Etapa 8 — SDF source abstraction: structured AND unstructured HEX8 input (Phase 3)

````text
ETAPA 8 — Pluggable SDF source: structured grid (trilinear) + unstructured HEX8 (FE shape functions)
(depends on Etapas 1-5, especially Etapa 2; INDEPENDENT of Etapas 6-7 — may be done before them)

GOAL: Support two kinds of input fields behind one narrow interface, so topology-optimization
results can come either from a structured grid (current path; trilinear interpolation, as in the
reference implementations) or from an unstructured conforming HEX8 FE mesh (isoparametric
shape-function interpolation). The mesher itself (lattice fill, warp, stencils, connectivity)
does NOT change — the generation lattice stays structured (that is a property of isosurface
stuffing); only the SDF SOURCE becomes pluggable.

WHY THIS WORKS (state these checks in the opening report):
  • The entire pipeline consumes the field through ONE function (eval_sdf) plus the cached
    node_sdf, so the seam is already narrow (verified 2026-06-10: TetGenerator :86/:103, Stencils
    centroid test :512; the only other access is get_cell_sdf_values for cell skip tests at
    TetGenerator :15/:46. NOTE: CorrectMeshVolume was deleted with volume correction, so it is no
    longer an eval_sdf call site).
  • HEX8 trilinear shape functions give a C0-continuous field on a CONFORMING mesh (the
    restriction to a shared face depends only on the 4 face nodes), so the zero isosurface is
    crack-free. Per-element interpolation without shared-face consistency would not be
    (the diagonal-choice problem) — hence the conforming-mesh requirement below.
  • After Etapa 2 the warp uses alpha = phi_i/(phi_i - phi_j) (scale-invariant) and no gradients,
    so the field only needs the right zero set + monotonicity near it — it does not have to be a
    true signed distance. SIMP density fields qualify.

IMPLEMENT:
  • New dir src/Fundamentals/SDFSources/ with:
      - SDFSource.jl: abstract type SDFSource; interface eval_sdf(src, p), bbox(src).
      - StructuredSDF.jl: wraps the current grid + values; MOVE the body of eval_sdf
        (trilinear interpolation, SDFOperations.jl:24) here unchanged.
      - UnstructuredSDF.jl: nodes, HEX8 connectivity, nodal phi. eval = point location
        (uniform spatial hash over element bounding boxes) -> inverse isoparametric mapping
        (3x3 Newton for (xi,eta,zeta), few iterations) -> trilinear shape functions.
        REQUIRE a conforming mesh (no hanging nodes) and positive Jacobians; validate on load.
      - Adapters.jl:
          SIMP: element-constant densities -> nodal values (volume-weighted average of adjacent
                elements), then phi = iso_level - rho with iso_level default 0.5, user-adjustable.
          Level-set: nodal field used directly; parameter for the input sign convention so the
                in-memory field is phi < 0 = inside (Etapa 1 convention).
          Auto-detect: given nodes + HEX8 connectivity, recognize a tensor-product grid (unique
                sorted x/y/z coordinates, complete lattice within tol) and build a StructuredSDF
                (fast path); otherwise UnstructuredSDF. Allow an explicit override parameter.
  • BlockMesh: add field sdf_source::SDFSource; keep ALL existing fields and the output contract
    (X, IEN, INE, node_sdf, node_hash) unchanged. mesh.SDF stays as the lattice-corner cache:
    for structured input it IS the input array (bit-identical behaviour); for unstructured input
    it is sampled from the source once at lattice construction.
  • Lattice decoupling: for structured input the lattice = the input grid (current behaviour,
    results must be IDENTICAL). For unstructured input build the lattice from bbox(source) +
    user-chosen dx + 2-cell padding (mirror quartet main.cpp's grid sizing).
  • eval_sdf(mesh, p) becomes a thin delegation to eval_sdf(mesh.sdf_source, p);
    get_cell_sdf_values keeps reading the lattice cache as today.
  • Outside-domain rule for unstructured sources: lattice points outside the hex mesh get a
    positive (outside) value (closest-element extrapolation or clamp). When the input is a
    design-domain box (the SIMP case), only the padding ring is affected.
  • Volume reference (CalcVolumeFromSDF — CorrectMeshVolume is gone): dispatch per source — current
    code for structured; Gauss quadrature over hex elements with FE interpolation for unstructured.
  • Clear explicit loops + English comments (CODE STYLE); the inverse mapping and the spatial
    hash deserve short explanatory docstrings.

VERIFY:
  1. Structured path: beam + gripper reproduce the Etapa-5 numbers EXACTLY (pure refactor there).
  2. Round-trip test: re-express the beam field as an unstructured HEX8 mesh (nodes+connectivity
     generated from the structured grid) and force it through the unstructured code path
     (auto-detect overridden) — the output mesh must match the structured run within
     floating-point tolerance. This validates point location + inverse mapping + shape functions.
  3. A genuinely unstructured small test: a distorted/graded hex block with an analytic field
     (e.g. sphere); check surface fidelity, watertightness, 0 inverted.
  NOTE on verification limits: quartet can only consume structured grids, so there is NO quartet
  oracle for unstructured inputs (the cross-validation section applies to paths 1-2 only). For
  unstructured cases judge by the quality histogram (the harness dih() metrics) + watertightness;
  once the owner approves a verified-good unstructured run, freeze its metrics as the automated
  regression baseline in test/.
CLOSEOUT (README: document the two input kinds, the adapters incl. iso_level, the conforming-mesh
requirement, and that the generation lattice itself remains structured by design).
````

---

## Etapa 9 — Functional-surface protection & optional AM-ready surface output (OPTIONAL; do LAST)

````text
ETAPA 9 — Protect functional surfaces from node-moving steps; optional AM surface smoothing
(OPTIONAL, lowest priority — a SCOPE EXPANSION beyond the quartet/Labelle volume-mesher goal.
Do only after Etapas 1-8, only if the use case calls for it. Independent of the core pipeline.)

MOTIVATION: Any step that MOVES surface nodes after meshing (volume correction, surface smoothing,
mesh optimization) risks displacing FUNCTIONAL surfaces — the regions where FE boundary conditions
live (loads, supports, symmetry planes, bolt holes, bearing faces). Moving those invalidates the
downstream FE/TO analysis. The mesher currently protects PLANAR functional regions via
plane_definitions + is_on_plane (manual, planar only). This stage generalizes that protection and,
optionally, adds a manufacturable (AM-ready) smoothed surface output.

KEY DECISIONS (from the 2026-06-10 assessment — keep in mind):
  • In a TO/FE workflow the functional surfaces are KNOWN UPSTREAM (they ARE the BC regions of the
    optimization problem). PREFER propagating that ground truth into the mesher (generalize
    plane_definitions to non-planar / tagged regions: cylinders for holes, arbitrary caller-supplied
    node/face tag sets) over heuristically re-detecting it from the faceted geometry.
  • Geometric AUTO-DETECTION earns its place ONLY when that info is lost/absent: third-party STL/SDF
    with no provenance, or when non-planar functional features must be found automatically.
  • Do NOT use functional-surface locking to justify keeping a flawed node-mover. (Volume correction
    is being removed for exactly its surface-moving / surface-fidelity problems.)

REFERENCE (external repo, not in Literature/):
  /Users/ondra/github/Disertace/Literature_md/2_State_of_the_Art/
      2021_Bacciaglia_SurfaceSmoothingForTopologicalOptimized.md
  An STL surface-smoothing method for TO→AM with two automatic "no-smoothing-space" detectors:
    - detect_flat_surface: cluster facets by shared normal (threshold L) -> planar regions;
    - detect_holes_edges: closed loops of ~90° sharp edges -> hole rims (perpendicular holes only).
  Both FREEZE the detected vertices during smoothing AND volume rescaling. NOTE: the smoothing CORE
  (HC-SDU) is NOT wanted for the FE-quality volume mesh (it fights the Labelle Hausdorff/dihedral
  guarantees); only the freeze-set concept transfers, unless an AM-ready surface output is added.

TASKS (pick per actual need; all OPTIONAL, default OFF):
  A. Generalize functional-region protection: extend plane_definitions / BoundedPlane to non-planar
     tagged regions, reusing boundary-face extraction. Any node-moving step consults this freeze-set.
  B. (Only if inputs lack BC provenance) Automatic functional-feature detection on the EXTRACTED
     boundary: robust flat-region detector (normal clustering with a tolerance suited to the faceted
     isosurface boundary) + hole-rim detector (sharp-edge loops). Validate against the known
     plane_definitions as ground truth on beam/gripper before trusting it.
  C. (Only if AM-ready output is in scope) An OPTIONAL feature-preserving surface smoothing pass
     (Bacciaglia HC-SDU or Taubin) on a COPY of the boundary surface, freeze-set held fixed, with
     volume rescaling — exported as a separate manufacturable STL, NEVER mutating the FE volume mesh.

VERIFY: feature OFF -> beam/gripper byte-identical to Etapa-8. Protection ON -> zero displacement on
the freeze-set. If smoothing is added: FE volume mesh unchanged, only the exported AM surface differs;
report compliance impact (the paper saw ~2% on the GE bracket). CLOSEOUT as usual (README: document
it as optional, default OFF, a scope expansion beyond the isosurface-stuffing core).
````
