# Performance TODO: stuffing with `cut_points = :bisection`

Found 2026-11-02 while timing stage 3 of the thesis recompute (gripper, session-5 clipped
field `min(f − c, d_Ω)`, lattice 49×86×128, ~2M tets, 1 thread; code of tag
`thesis-recompute-s5`). Not fixed yet.

## Measurement

| stage | time [s] | `eval_sdf` calls |
|---|---:|---:|
| lattice corners (`BlockMesh`) | 0.9 | 0.54 M |
| A15 fill (`generate_mesh!`) | 19.3 | 12.7 M |
| **`warp!` (bisection)** | **210.9** | **55.9 M** |
| slicing (bisection) | 10.5 | 2.6 M |
| connectivity + inverted + isolated | 7.0 | 0 |
| plane alignment (`warp_mesh_by_planes_sdf!`) | 14.2 | 0 |
| **total** | **262.7** | 71.7 M |

The same run with the plain RBF field (no clip) takes 39.6 s: the call counts are the same,
one call costs 0.28 µs instead of 3.8 µs near the surface (the `HexDomain` part of the clip
dominates). The stuffing itself is cheap; the time is redundant field evaluation.

## Causes

1. `warp!` bisects every crossing edge once per incident tet
   (`src/GenerateMesh/TetGenerator.jl`, the loop over `mesh.IEN`): 1,397,564 visits of
   274,644 unique crossing edges (×5.09), 40 `eval_sdf` each → 56 M calls instead of 11 M.
   The comment "keeping the per-vertex minimum makes the repetition harmless" is only true
   for the linear `alpha` (quartet); the docstring of `bisect_cut_alpha` ("40 calls per
   crossing edge") undercounts it by the same factor.
2. `slice_ambiguous_tetrahedra!` / `cut_edge!` bisect again edges that `warp!` already
   bisected (their endpoints were not moved, so `alpha` is the same): 2.6 M calls, 10 s.
3. The C-F1 cell-skip test in `process_cell_A15!` evaluates all 27 tile points of every
   cell whose corners and face neighbours are outside, without a cache: 447,325 cells (86 %)
   × 27 = 12.08 M calls for 3.72 M unique points (×3.25). On the gripper it keeps no cell.
4. `find_surface_nodes` (`src/Modification/ModifyResultingMesh.jl`) keys a `Dict` by
   `Set{Int}` (8 M `Set` allocations on 2M tets): 8–9 s of the 14 s alignment.

## Reference: quartet (`make_tet_mesh.cpp`)

The A15 tile spans 4×4×4 SDF grid cells, so every lattice vertex is a grid node: `phi` is
read from the stored grid once per vertex (`lattice_map`), the keep test runs per tile tet
on stored values, and `warp_vertices` / `cut_edge` use the linear
`alpha = phi_i / (phi_i − phi_j)`. Nothing is evaluated after the lattice is built, which is
why stuffing is cheap there. Here the tile spans one grid cell (nodes at h/4) and the
`:bisection` mode evaluates the field along edges, so the evaluation count matters.

## Plan (estimates for the gripper case above)

| change | effect | output |
|---|---|---|
| cache `alpha` per canonical edge in `warp!` | 56 M → 11 M calls, ~211 s → ~41 s | identical up to last bits (edge orientation) |
| pass that cache on to the slicing | −10 s | identical |
| C-F1: evaluate each A15 lattice vertex once (keyed by its integer lattice coordinate, like `mesh.node_hash`), reuse it for `node_sdf` at node creation | 12 M → 3.7 M calls, −13 s | identical |
| `find_surface_nodes` with `face_key` tuples | −8 s | identical |
| optional: Brent/Illinois instead of 40 halvings (smooth field, ~6–8 calls) | `warp!` → ~8 s | same root to ~1e-12 |

Without changing the output: ~263 s → ~60 s; with the faster root finder ~25–30 s.
Profiling scripts used for the numbers: counting `SDFSource` wrapper around the stage-3
source, stage-by-stage `@elapsed` (thesis recompute session, 2026-11-02).
