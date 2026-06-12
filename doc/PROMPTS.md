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
  • Etapas 1-5 — DONE, committed on `dev` (Etapa 4 = 25c6f70, Etapa 5 = 0bc8409). Conformance to the
    published quartet/Labelle algorithm verified (audit Part A); considered reviewed.
  • Post-audit performance cleanups (audit Part C) — DONE: slice/trim allocation cuts (95871b7),
    candidate-only face map C1+C2 (96804f8), build INE once C3 (ce38be4). Gripper ≈14.5s → ≈6.0s.
  • Test suite overhaul — DONE (a13232f): `Pkg.test()` is the real gate now (64/64) — invariants +
    frozen baselines + per-stage beam VTU diagnostics; examples moved to top-level examples/.
  • cut_points = :bisection — DONE (0af0d3b): optional L&S §3.1 edge-bisection mode for non-distance
    (smoothed) SDF inputs; default stays :linear (quartet-faithful).
  • Volume correction — REMOVED (c224f05); the pipeline has no correct_mesh_volume! step, so any
    "compare before correct_mesh_volume!" note below is moot.
  • Etapa 9 — only the prompt text lives in this file (12e9bc9); no code yet.
  • Etapa 10 (gated mesh relaxation: :uniform size-equalizing springs + :quality quartet-style
    maximin smoothing, one switch, fixed topology, default OFF) — DONE 2026-06-12 (uncommitted at
    time of writing; src/Modification/RelaxMesh.jl, wired into MeshGenerationOptions.relax). DESIGN
    DEVIATION from this prompt, approved by results: the gate metric is a TWO-SIDED check on the
    dihedral-angle cosines (sharpest angle may not sharpen AND widest may not widen), NOT the eta
    mean-ratio the prompt suggested for :uniform — eta does not bound the dihedral, so it would not
    make "min dihedral >= pre-pass" a true construction; the two-sided cosine gate makes BOTH "min
    does not drop" and "max does not grow" hard guarantees for both modes. Results (default
    :uniform): beam min 11.28->18.37, max 153.21->131.62, >140 81->0; gripper min 9.54->12.35,
    <10 11->0, >140 1136->91. OFF path byte-identical (2392/2392 default, 2416/2416 +gripper).
    Determinism verified. Timing ~1s beam, ~12-15s gripper (thin walls -> ~2/3 of nodes in band).
    Independent of Etapas 6-7/9; needs 1-5 + 8.
  • CURRENT FOCUS: **Stabilization / review** before Phase 2 — close the audit's Part B (dead code,
    never executed; only Part C was). Confirmed candidates: orphan src/Fundamentals/Hex8_shape.jl
    (not include-d, shape_functions unreferenced); compute_gradient (Newton-warp leftover, no caller,
    still exported in Fundamentals.jl); over-structured CaseParams (only nnzz dihedral bounds are read,
    the nzzz field is dead). DELETE NOTHING without owner sign-off (audit rule).
  • NEXT coding stage (Phase 2, after stabilization): **Etapa 6** (optional adaptive octree sizing).
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
  src/GenerateMesh/DihedralAngles.jl    — DihedralBounds (:17), create_warping_params (:34), compute_dihedral_angle_range (:50)
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

## 🔎 Current-state audit — conformance / dead code / performance (READ-ONLY review)

Not a numbered stage and NOT a coding task. This prompt produces a written REPORT on the state of
the code after Etapas 1-5, as a gate before Etapa 6 is started. It changes NO code: it lists
findings with evidence and recommendations, then stops for owner review. Approved findings become
their own scoped cleanup stages afterwards (with the usual closeout). Paste the **Common context**
block first, then this prompt.

````text
CURRENT-STATE AUDIT — Implicit2TetMesh.jl (READ-ONLY: produce a report, change no code, do not commit)

GOAL: Assess the post-Etapa-5 codebase along three axes — (A) conformance to the published/quartet
algorithm, (B) dead code left over from the pre-Etapa implementation, (C) performance — and report
findings with EVIDENCE (file:line, measured numbers) and prioritized, low-risk recommendations.
Do NOT edit code, do NOT commit. End by STOPPING for review.

SCOPE NOTE — what "conform to the original" means here: the TARGET is the PUBLISHED algorithm, i.e.
quartet's warp_vertices + trim_spikes AND Labelle §3.4 for surface tets. Etapa 4 DELIBERATELY
replaced quartet's centroid-only remove_exterior_tets with Labelle §3.4 (more principled), so the
FINAL mesh is intentionally NOT bit-identical to quartet's final output. Do not flag that documented
divergence as a bug. The warp and the trim stencils, by contrast, ARE meant to match quartet 1:1.

BEFORE STARTING: follow the common-context "BEFORE STARTING" checklist (read the listed code +
references). Re-grep every function for its CURRENT line number — the "Key files" numbers in this
file have drifted. Read quartet src/make_tet_mesh.cpp: warp_vertices (~152-193), trim_spikes
(~221-346), remove_exterior_tets (~350-366), and predicates.cpp (sign convention). Skim the project
auto-memory notes (quartet-crossval-harness, watertightness-pinches, volume-correction-removed) for
the verified Etapa-2..5 results and the orient() sign GOTCHA.

PART A — CONFORMANCE TO THE PUBLISHED / QUARTET ALGORITHM
  A1. warp! (TetGenerator.jl) vs quartet warp_vertices: edge-based, LINEAR cut point (no Newton),
      alpha = phi_i/(phi_i - phi_j), threshold 0.3, each vertex keeps its CLOSEST qualifying cut
      point, displacements applied at once (order-independent), warped node_sdf set to 0. Confirm 1:1.
  A2. slice_ambiguous_tetrahedra! / apply_stencil_trim_spikes! (Stencils.jl) vs trim_spikes: every
      sign case (+++-/+---/++--/+0--/++0-/+00-), the QUAD-SPLIT diagonal choice (the sorted-pqrs /
      vertex-index tie-break that keeps face-adjacent tets watertight), and orientation handling.
  A3. Surface-only ("quadruple-zero") handling vs Labelle §3.4: inverted-or-bad-dihedral → discard;
      then 4 faces adjoin solid → retain / 0 adjoin → discard / else centroid SDF sign. Verify it
      matches the PUBLISHED heuristic (this is the INTENDED divergence from quartet — judge it correct,
      not whether it equals quartet).
  A4. Sign convention: phi<0=inside applied consistently (Etapa 1); exact predicates back EVERY
      orientation/inversion decision through the single is_positively_oriented helper, with the
      "orient() returns the OPPOSITE sign of the float det" convention correct (memory GOTCHA). A
      tolerance gates ONLY near-zero-volume magnitude, never the sign.
  A5. A15-only: no Schläfli, no experimental_nzzz remnants in any ACTIVE code path.
  A6. RE-RUN the quartet cross-validation (🧪 section above) on beam at upsample factor 4 (the only
      factor where the two lattices physically coincide): rebuild the driver, export the SDF, run,
      compare. Expect the post-warp mesh still bit-identical and the final histograms as recorded in
      memory. If the C++ build is unavailable, FALL BACK to scratch/measure.jl vs the frozen baselines
      (beam: min 11.276° / <10°=0 / >140°=81 / inverted 0; gripper: min 9.537° / <10°=11 / >140°=1136
      / inverted 0) and state that the oracle run was skipped.
  Report every divergence from the published algorithm as a finding (severity + file:line + fix idea).

PART B — DEAD CODE / LEFTOVERS OF THE PRE-ETAPA IMPLEMENTATION
  Method: build the reachable call graph from the TWO real entry points — generate_tetrahedral_mesh
  (TetMeshGenerator.jl) and the scripts under test/ (runtests.jl, Examples/beam.jl, Examples/gripper.jl).
  A symbol reachable from neither, and not part of the package's intended public API, is a candidate.
  Check, at minimum:
   • Files in src/ that are never include-d (suspect: Fundamentals/Hex8_shape.jl — it is NOT listed in
     Fundamentals.jl; confirm). List every such orphan file.
   • NewCases-Experimental.jl: only create_warping_params + compute_dihedral_angle_range are meant to
     survive — list every other definition there and whether anything references it.
   • Confirmed-removed features that might leave danglers: Newton warp (warp_node_to_isocontour! — and
     is compute_gradient now called by anything?), Schläfli (process_cell_Schlafli!, SchlafliScheme.jl),
     post-hoc volume correction (correct_mesh_volume! / CorrectMeshVolume / unused parts of
     CalcVolumeFromSDF). Confirm each is fully gone, or flag the remnant.
   • Utils/ and Modification/ exports — assess_mesh_quality, slice_mesh_with_plane!,
     count_negative_determinants, calculate_volume_from_sdf, ModifyResultingMesh internals: which are
     reachable from the pipeline OR from tests/examples? Test-only or public-API use COUNTS as used — say so.
   • BlockMesh struct fields never read (e.g. node_hash, grid_tol?) and MeshGenerationOptions
     fields/params that no longer change the output.
   • module export lists naming symbols that nothing imports.
  For EACH candidate, classify: (1) safe to delete, (2) keep — public API or test/example only,
  (3) uncertain — needs an owner decision. Produce a removal CHECKLIST. DELETE NOTHING in this pass.

PART C — PERFORMANCE
  Use the existing harnesses (scratch/, see memory): timeit.jl (per-phase timing, JIT-warmed on beam),
  slice_probe.jl (alloc/GC per phase), measure.jl (quality). Re-measure gripper end-to-end + per phase.
   • Compare to the Etapa-5 baseline (gripper ≈8.1s total). Confirm NO regression; give a per-phase
     timing + allocation table.
   • Validate (or refute, with data) the recorded claim that slice's remaining ≈2.5s is GENUINE
     trimming work (cut-point creation, output element Vectors, exact orient per output tet), not waste.
   • Confirm the 3 update_connectivity! calls are each still REQUIRED (memory: all three are).
   • Look ONLY for NEW, low-risk, OUTPUT-PRESERVING opportunities (allocation hotspots, redundant
     eval_sdf, avoidable rebuilds). For each: expected impact + risk. Do not re-litigate solved wins;
     do not propose algorithm changes for speed in this audit.

DELIVERABLE: ONE written report, three sections (A/B/C). Each finding: what it is, evidence (file:line
or measured number), severity/priority, and a concrete recommendation. No code changes, no commit.
STOP for review. (If asked afterwards, the cleanups become their own scoped stages with the usual closeout.)
````

---

## 🧰 Test suite overhaul — functional tests + per-stage beam diagnostics

Not a numbered stage: an infrastructure task that can run any time after Etapa 5. It replaces the
current script-style `test/runtests.jl` (zero `@test` assertions) with a real, assertion-based test
suite, plus a per-stage beam diagnostic that exports a VTU after every pipeline phase so a broken
phase can be SEEN in ParaView. Paste the **Common context** block first, then this prompt.

````text
TEST SUITE OVERHAUL — Implicit2TetMesh.jl (replace script-style runtests.jl with real functional tests)

GOAL: Make `Pkg.test()` a meaningful gate. Two kinds of tests:
  (a) INVARIANT tests — properties EVERY correct output must satisfy (watertight boundary, zero
      inverted tets, single component, node-SDF consistency, positive volumes). These survive future
      Etapas 6-9 unchanged.
  (b) BASELINE (regression) tests — exact frozen numbers (node/tet counts, dihedral histogram). The
      pipeline is deterministic (byte-identical reruns verified in Etapa 5), so exact equality is
      safe and sharp. Baselines live in ONE file so an INTENTIONAL pipeline change (e.g. Etapa 6)
      re-freezes them in a single, reviewable commit.
Plus a per-stage beam diagnostic test that exports the mesh to VTU after EACH pipeline phase, so
when a stage breaks, the assertion names the phase and the VTU shows the damage.

TARGET LAYOUT (owner requirement — keep it clean; the TRACKED tree under test/ is EXACTLY two
files and one folder):
  test/
    runtests.jl            — the core suite itself (all @testsets for (a)+(b) live HERE, not in
                             per-test include files), plus `include("test_beam_stages.jl")`
    test_beam_stages.jl    — the per-stage beam diagnostic (also runnable standalone)
    helpers/               — ALL helper code, no @testsets inside:
      helpers.jl           — check helpers (see below) + data-loading convenience
      baselines.jl         — the frozen constants, nothing else
  test/output/             — runtime VTU artifacts ONLY; created by the tests at runtime, added to
                             .gitignore (it is not part of the tracked tree, which keeps the
                             two-files-one-folder rule for what is in git)

RESTRUCTURING (to reach that layout):
  1. Move `test/Examples/beam.jl` and `test/Examples/gripper.jl` to a new top-level `examples/`
     directory (they are documentation, not tests — standard Julia convention). Fix their relative
     data paths for the new location; delete `test/Examples/`.
  2. Fold `test/GenerateMeshTests/validate_sdf_values.jl` into `test/helpers/helpers.jl`: strip the
     println banners (make it silent), have it RETURN the stats Dict only; delete
     `test/GenerateMeshTests/`.
  3. Delete the stray `.vtu` files lying in `test/` (untracked build artifacts). All test exports go
     to `test/output/` from now on; add `test/output/` to .gitignore.
  4. Use `joinpath(@__DIR__, "..", "data", ...)` for all data paths — never pwd-relative strings.

HELPERS (test/helpers/helpers.jl — promote the proven logic from scratch/, do not reinvent):
  • check_watertight(mesh) -> (open_edges, boundary_faces, ...): port from scratch/watertight.jl.
    CORRECTNESS CRITERION IS open_edges == 0 ONLY. Do NOT assert non-manifold pinch counts — they
    are benign thin-feature geometry and NOT invariant across stages (see memory: beam 1→7 after
    Etapa 4, expected).
  • dihedral_stats(mesh) -> (min, max, <5°, <10°, >140°, inverted): port from scratch/measure.jl
    (the dih() function), unchanged semantics so numbers stay comparable with all recorded history.
  • count_inverted_exact(mesh): use the SAME exact-predicate orientation helper the pipeline uses
    (is_positively_oriented). GOTCHA (memory): ExactPredicates.orient returns the OPPOSITE sign of
    the float det — reuse the pipeline helper, never call orient() raw.
  • count_components(mesh): connected components over shared faces (reuse the face->elements map
    pattern from RemoveIsolatedComponents.jl).
  • validate_node_sdf_values(mesh, tol): the silenced version from step 2 above.
  • load_beam(), load_gripper(): JLD2 loading with the @__DIR__-based paths.

CORE SUITE (test/runtests.jl — clean nested @testsets, no RUN_* flags, no println noise):
  @testset "Options validation": MeshGenerationOptions rejects scheme != "A15" and negative
    warp_param (@test_throws); defaults are as documented.
  @testset "Beam pipeline (no planes)": run generate_tetrahedral_mesh on beam, output prefix under
    test/output/. Assert ALL invariants: open_edges == 0, inverted == 0, components == 1,
    validate_node_sdf_values max_error <= 0.005, all element volumes > 0. Then assert BASELINES:
    exact node count, exact tet count, dihedral histogram fields == frozen constants.
  @testset "Beam pipeline (with cutting planes)": the two-plane beam config (Square(30) at x=0,
    Square(5) at x=60, warp_param 0.3). After the cut: open_edges == 0, inverted == 0, and the
    baseline counts for the cut mesh.
  @testset "Beam per-stage diagnostics": include("test_beam_stages.jl").
  @testset "Gripper (opt-in)": ONLY when ENV["I2TM_TEST_GRIPPER"] == "1" — too slow for the default
    run. Same invariants + frozen gripper baselines. Default suite must stay fast (~1-2 min).

PER-STAGE BEAM DIAGNOSTIC (test/test_beam_stages.jl):
  Purpose: when a pipeline phase breaks, the failing assertion NAMES the phase and the exported VTU
  SHOWS it. Run the beam pipeline stage by stage with the EXACT call sequence of
  generate_tetrahedral_mesh (src/TetMeshGenerator.jl — keep the two update_connectivity! calls with
  build_ine = false; if the sequence there ever changes, mirror it here).
  After EVERY stage: FIRST export the mesh (export_mesh_vtu) to test/output/ with a numbered name —
  Beam-stage_1-generated.vtu, _2-warped, _3-sliced, _4-inverted_removed, _5-components,
  _6-plane_cut — THEN run that stage's assertions, in this order, so the VTU exists for inspection
  even when the assertion fails. One @testset per stage:
    Stage 1 generate_mesh!: nodes > 0, tets > 0, every IEN index in bounds, no NaN/Inf coordinates,
            count_inverted_exact == 0.
    Stage 2 warp!: no NaN coordinates, validate_node_sdf_values max_error <= 0.005 (warped nodes
            carry node_sdf = 0 and must sit on the isosurface within tolerance).
    Stage 3 slice_ambiguous_tetrahedra! (+ update_connectivity!): open_edges == 0 (the watertight
            check), every IEN index in bounds.
    Stage 4 remove_inverted_elements!: count_inverted_exact == 0.
    Stage 5 remove_isolated_components!(keep_largest) (+ final update_connectivity!): components
            == 1, open_edges == 0, dihedral histogram == frozen stage-5 baseline (this stage equals
            the no-planes pipeline output, so the constants are shared with the core suite).
    Stage 6 warp_mesh_by_planes_sdf! (beam two-plane config): open_edges == 0,
            count_inverted_exact == 0.
  The file must also work standalone (julia --project=. test/test_beam_stages.jl) for interactive
  debugging — guard the include-vs-standalone difference cleanly (e.g. it defines/uses the helpers
  via include of helpers/ files that are idempotent to include twice).

BASELINES (test/helpers/baselines.jl):
  Re-freeze ALL constants from the CURRENT dev HEAD at implementation time by running the pipeline
  once and copying the numbers — do NOT trust the numbers in this prompt to still be current. For
  reference, the recorded post-Etapa-4/5 values were: beam tets 92,532, min 11.276°, <10° 0,
  >140° 81, inverted 0; gripper tets 1,984,316, min 9.537°, <10° 11, >140° 1136, inverted 0.
  Each constant gets a one-line comment. Top of file: a short note that after an INTENTIONAL
  pipeline change these are re-frozen by running the suite and updating this one file in the same
  commit as the change.

OUT OF SCOPE (do not do): quartet C++ cross-validation in tests (stays in scratch/); soft quality
thresholds (min angle > X°) — only invariants and exact baselines; synthetic/analytic SDF fixtures
(wait for Etapa 8); touching src/ beyond what the moved helpers need (they must need NOTHING —
if a src change seems required, STOP and report); Project.toml [extras]/[targets] reshuffling.

VERIFY: `julia --project=. -e 'using Pkg; Pkg.test()'` passes (fast path). Then once with
I2TM_TEST_GRIPPER=1 — passes. Break one stage on purpose locally (e.g. skip warp!) and confirm the
per-stage test fails AT THE RIGHT STAGE and the VTUs up to that stage exist in test/output/. Revert.
Confirm `git status` shows a clean test/ tree: exactly runtests.jl, test_beam_stages.jl, helpers/.

CLOSEOUT: as in the common context — but the regression harness step is now `Pkg.test()` itself
(scratch/measure.jl remains for ad-hoc gripper measurements). Update README (testing section: how to
run, the gripper ENV flag, where VTUs land, examples/ move). Commit to dev. STOP for review.
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

---

## Etapa 10 — Gated mesh relaxation & quality optimization (OPTIONAL post-pass; fixed topology)  ✅ DONE (awaiting review)

````text
ETAPA 10 — Optional gated vertex relaxation: mode :uniform (size equalization, "springs") and
mode :quality (quartet-style maximin smoothing), behind one shared skeleton and one user switch.
(OPTIONAL, default OFF. Depends on Etapas 1-5 + 8 (final INE build, eval_sdf seam, exact
predicates, pluggable sources); INDEPENDENT of Etapas 6-7 and 9 — may be done before them.
FORWARD-COMPAT: this pass MOVES nodes, so when Etapa 9's functional-surface freeze-sets land,
relaxation MUST consult them; design the freeze-set as an input parameter from day one.)

MOTIVATION: the A15 interior is already uniform and near-optimal by construction; ALL element-size
variance and ALL bad dihedral angles live in the 1-2 element layers created by warp + trim. Two
different (sometimes conflicting) post-processing objectives follow:
  • :uniform — equalize edge lengths near the boundary (FE users want similar-sized elements);
  • :quality — actively lift the worst dihedral angles (gripper today: min 9.537°, 11 tets <10°,
    1136 >140° — those tails are the target).
Both fit one skeleton: PROPOSE a new vertex position -> GATE it on local element quality ->
PROJECT surface vertices back to phi = 0. They differ only in the proposal generator.

DECISIONS ALREADY MADE (do not re-litigate):
  • FIXED TOPOLOGY ONLY: vertex relocation, never flips/collapses/insertions (Stellar / CGAL
    tetrahedral remeshing are out of scope — they would break the watertight trim structure,
    determinism, and the runtime budget).
  • DEFAULT OFF. With relaxation off, beam/gripper output stays byte-identical; the quartet
    cross-validation oracle is only ever compared with relaxation off.
  • Surface vertices STAY ON phi = 0 at all times: tangential proposal + re-projection by
    BISECTION ALONG THE VERTEX PSEUDO-NORMAL. NEVER project along the SDF gradient (quartet does;
    it only works for true distance fields — ours may be SIMP/level-set; compute_gradient was
    deleted in audit Part B, do NOT re-add it). This is NOT a volume corrector: it never moves
    nodes off the zero level (the removed correct_mesh_volume! failure mode), so the README
    volume-accuracy contract stays intact. Volume changes only by chord redistribution.
  • DETERMINISTIC: no rand() (quartet's optimizer uses rand — do not port that). Fixed-order
    Gauss-Seidel sweeps (or Jacobi), fixed candidate patterns. Same input -> same output, so
    relaxed baselines can be frozen in test/helpers/baselines.jl.
  • The gate makes safety a CONSTRUCTION, not a hope: the current position is always a candidate,
    an accepted move never lowers the local min quality => global min quality is monotonically
    non-decreasing => "min dihedral >= pre-pass" is a HARD testable assertion.

BEFORE STARTING: follow the common-context checklist. Read:
  • Literature/quartet-original/src/optimize_tet_mesh.cpp (the whole file, ~265 live lines) — the
    in-family precedent: per-vertex candidate set incl. the original point, maximin over incident
    tets, tie -> smallest move, interior/boundary/feature vertex classes, early stop;
    tet_quality.cpp (min-sine-of-dihedral metric, orientation-signed) and sdf.cpp
    projectToIsosurface (bracket clamped to +-dx — keep the clamping idea, replace gradient with
    pseudo-normal).
  • Published basis to confirm and cite in the opening report: Persson & Strang 2004 "A Simple
    Mesh Generator in MATLAB" (DistMesh spring force = the :uniform proposal); Persson 2006 "Mesh
    size functions for implicit geometries and PDE-based gradient limiting" (curvature sizing +
    |grad h| <= g); Freitag 1997 (smart Laplacian = propose-then-gate); Thürmer & Wüthrich 1998 /
    Bærentzen & Aanæs 2005 (angle-weighted pseudo-normals).
  • src/Modification/ModifyResultingMesh.jl find_surface_nodes / surface_faces (count==1 boundary
    extraction to reuse) and Stencils.jl is_positively_oriented (the ONLY orientation check to
    use; remember ExactPredicates.orient's sign is OPPOSITE the float det).

IMPLEMENT (new file src/Modification/RelaxMesh.jl, exported relax_mesh!(mesh, opts); explicit
loops + English docstrings per CODE STYLE):
  • RelaxOptions struct (all defaults are starting points — tune while measuring):
      mode::Symbol = :uniform | :quality
      max_sweeps::Int (10), tol: stop when max displacement < 1e-2*grid_step (:uniform) or when
        min quality stops improving (:quality, quartet's early-stop idea)
      band::Int = 2 — active set = surface vertices + `band` topological rings inward (the A15
        interior is already optimal; band = typemax(Int)/:all as an override for experiments)
      omega::Float64 = 0.3 (:uniform step factor)
      sizing::Symbol = :uniform | :curvature (only meaningful for mode = :uniform), with
        curvature params: alpha, h_min/h_max clamps, grading g = 0.3
      frozen::AbstractVector{Int} = [] — externally supplied freeze-set (Etapa 9 hook)
    Wire into MeshGenerationOptions as relax::Union{RelaxOptions,Nothing} = nothing (OFF).
  • Shared infrastructure (built once per call):
      - boundary triangulation: faces with incidence 1 (reuse the count==1 pattern);
      - angle-weighted pseudo-normals at surface vertices from that triangulation;
      - FREEZE: vertices on non-manifold boundary edges (pinches — thin-feature sheets, see
        watertightness lore: their pseudo-normal is meaningless), vertices in opts.frozen, and
        anything outside the active band;
      - vertex->neighbor lists from IEN restricted to the band (INE already exists — the pass
        runs after the final update_connectivity!, and since topology never changes, INE, IEN
        and the boundary triangulation stay VALID throughout; only positions move).
  • The GATE (shared by both modes): a candidate position x' for vertex v is accepted only if
    over all tets incident to v (via INE): (a) is_positively_oriented holds for every tet (exact
    sign, no tolerance), (b) min quality does not decrease vs the current position, and (c) for
    INTERIOR vertices eval_sdf(mesh, x') < 0 (a relaxed interior node must stay strictly inside;
    keeps the node-SDF invariant meaningful and protects thin features). On failure try step
    s in {1, 1/2, 1/4} (:uniform), else keep the current position.
  • mode :uniform (the springs): proposal d_v = omega * sum_j (x_j - x_v) * (1 - L0/|e_vj|)
    over neighbor edges — compression for short edges, tension for long ones; equilibrium =
    uniform lengths (DistMesh force). L0 = median active-band edge length (sizing = :uniform).
    SURFACE vertices: neighbors = surface neighbors only (boundary triangulation), project the
    proposal into the tangent plane (d -= (d.n)n with n = pseudo-normal), apply, then re-project
    onto phi = 0 by bisection of eval_sdf along the pseudo-normal, bracket clamped to
    +-0.5 * local edge length (avoid jumping to the opposite sheet of a thin wall); the vertex's
    node_sdf stays exactly 0.0 by construction. Gate metric: mean ratio
    eta = 12 * (3V)^(2/3) / sum(l^2) (cheap, no trig; 1.0 = regular tet), sign from (a).
  • mode :quality (quartet, deterministically): candidate set per vertex = current position +
    a FIXED pattern — interior: 8 cube corners + 6 axis points at radius P = 0.5*grid_step, plus
    the same 14 at P/2; surface: the analogous 2D pattern (4 corners + 4 axes at P and P/2) in
    the tangent plane, each candidate re-projected to phi = 0 via the SAME bisection helper.
    Score = min over incident tets of MIN-SINE-OF-DIHEDRAL (quartet's metric — it is exactly the
    quantity our histograms report, so improvements are directly visible); pick the argmax,
    tie -> smallest displacement (quartet's tie-break). The gate is built in (current position
    competes). sizing/L0 plays no role in this mode.
  • Pipeline placement (TetMeshGenerator.jl): after remove_isolated_components! + the final
    update_connectivity! (INE just built), BEFORE TetMesh_volumes/export and BEFORE the optional
    warp_mesh_by_planes_sdf! (plane alignment must run last so it is not un-done; until Etapa 9
    lands, plane-region nodes are simply re-aligned afterwards as today).
  • Bookkeeping (easy to forget):
      - refresh node_sdf for every MOVED interior vertex (node_sdf[v] = eval_sdf(mesh, X[v]));
        surface vertices keep 0.0;
      - mesh.node_hash / node_map are coordinate-keyed merge artifacts — verify nothing consumes
        them after this stage (they are connectivity-merge tools); document, or empty! them;
      - count accepted/rejected moves and report per sweep (like warp!'s "Warped N vertices").
  • sizing = :curvature (the opt-in concentration feature, :uniform mode only): per-surface-vertex
    target L0(x) = clamp(alpha / kappa(x), h_min, h_max), kappa estimated DISCRETELY as the max
    angle between the vertex pseudo-normal and its surface-neighbors' pseudo-normals divided by
    edge length (no SDF derivatives); then GRADIENT LIMITING: one Dijkstra-like relaxation pass
    over the active-band edge graph enforcing h_j <= h_i + g*|e_ij| (discrete |grad h| <= g,
    Persson 2006 — this is the published answer to "denser spots must not pull and degrade their
    neighbors"); interior L0 = limited h propagated inward. Honest expectation (document it):
    with fixed topology this yields MILD concentration (edge ratios ~1.5-2x), not true
    refinement — real refinement is Etapa 6's octree.
  • OUT OF SCOPE (explicit follow-ups, do not build now): mode = :both (a :uniform pass followed
    by :quality sweeps — safe to chain since the gate is monotone, but ship the two modes first);
    worklist/priority-queue scheduling (only re-visit vertices whose star changed) — add only if
    profiling shows the plain band sweeps are too slow.

VERIFY (run scratch/timeit.jl before/after; budget: :uniform ~1 s, :quality a few s on gripper):
  1. OFF: full test suite green unchanged (64/64 today), beam + gripper byte-identical to the
     pre-Etapa-10 head. The quartet oracle workflow is untouched.
  2. :uniform ON (beam, default band/sweeps): watertight (open edges 0), 0 inverted (exact
     predicate), single component, every surface vertex has |eval_sdf| <= bisection tol, every
     interior vertex eval_sdf < 0, min dihedral >= the pre-pass value (sharp gate assertion),
     and the coefficient of variation of active-band edge lengths strictly DECREASES (this is
     the point of the mode — assert it).
  3. :quality ON (beam): same invariants, plus min dihedral strictly INCREASES vs pre-pass and
     the >140 deg count does not grow (expect it to shrink). Freeze the achieved numbers as the
     relaxed regression baseline AFTER the owner approves the first verified run.
  4. Determinism: run each mode twice — identical output (hash X and IEN).
  5. Unstructured source: the analytic-sphere HEX8 case (test_sdf_sources.jl) with relaxation ON
     must pass the same invariants; judge surface fidelity by GEOMETRIC distance to the true
     sphere — boundary_max_abs_sdf is NOT a valid metric for unstructured sources (known gotcha).
  6. Curvature sizing ON: edge lengths near high-curvature regions shrink toward alpha/kappa,
     the h-field respects the grading bound g, and ALL invariants of (2) still hold.
  7. With plane_definitions set: planes still come out aligned (plane warp runs after relax).
  8. OPTIONAL qualitative cross-check: build the quartet driver with optimize=true and compare
     the achieved min dihedral on beam f4 against our :quality mode. NOT bit-comparable (quartet
     uses rand(), float32 .tet output, gradient projection) — compare the quality statistics only.
CLOSEOUT as usual. README: document RelaxOptions (both modes, default OFF, the curvature option +
its honest limits), state explicitly that surface nodes remain on phi = 0 (why this is NOT the
removed volume correction), and update the TODO list. Update this file's Status block.
````
