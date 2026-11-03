# Performance: stuffing with `cut_points = :bisection`

Found 2026-11-02 while timing stage 3 of the thesis recompute (gripper, session-5 clipped
field `min(f − c, d_Ω)`, lattice 49×86×128, ~2M tets, 1 thread; code of tag
`thesis-recompute-s5`). Items 1–3 below are fixed (2026-11-03) without changing the output;
item 4 and the faster root finder are open.

## Measurement before the fix

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

1. `warp!` bisected every crossing edge once per incident tet
   (`src/GenerateMesh/TetGenerator.jl`, the loop over `mesh.IEN`): 1,397,564 visits of
   274,644 unique crossing edges (×5.09), 40 `eval_sdf` each → 56 M calls instead of 11 M.
   The comment "keeping the per-vertex minimum makes the repetition harmless" was only true
   for the linear `alpha` (quartet); the docstring of `bisect_cut_alpha` ("40 calls per
   crossing edge") undercounted it by the same factor.
2. `slice_ambiguous_tetrahedra!` / `cut_edge!` bisected again edges that `warp!` already
   bisected (their endpoints were not moved, so `alpha` is the same): 2.6 M calls, 10 s.
3. The C-F1 cell-skip test in `process_cell_A15!` evaluated all 27 tile points of every
   cell whose corners and face neighbours are outside, without a cache: 447,325 cells (86 %)
   × 27 = 12.08 M calls for 3.72 M unique points (×3.25). On the gripper it keeps no cell.
4. `find_surface_nodes` (`src/Modification/ModifyResultingMesh.jl`) keys a `Dict` by
   `Set{Int}` (8 M `Set` allocations on 2M tets): 8–9 s of the alignment.

## Reference: quartet (`make_tet_mesh.cpp`)

The A15 tile spans 4×4×4 SDF grid cells, so every lattice vertex is a grid node: `phi` is
read from the stored grid once per vertex (`lattice_map`), the keep test runs per tile tet
on stored values, and `warp_vertices` / `cut_edge` use the linear
`alpha = phi_i / (phi_i − phi_j)`. Nothing is evaluated after the lattice is built, which is
why stuffing is cheap there. Here the tile spans one grid cell (nodes at h/4) and the
`:bisection` mode evaluates the field along edges, so the evaluation count matters.

## Fixed (2026-11-03): items 1–3, output bit-identical

- **Items 1 and 2.** `cached_bisect_cut_alpha!` keeps every bisection result in a
  `CutAlphaCache`, keyed by the two end points of the oriented edge and `phi_i`, i.e. by
  everything `bisect_cut_alpha` reads. `generate_tetrahedral_mesh` passes one cache to
  `warp!` and to `slice_ambiguous_tetrahedra!`. The key holds no node indices, so the cache
  survives the renumbering of `update_connectivity!` between the two stages.
- **Item 3.** `tile_point_sdf!` keeps the SDF of every A15 tile point in a `TileSDFCache`,
  keyed by its integer lattice coordinate (the key of `mesh.node_hash`). The C-F1 test and
  the node creation read it. A stored value is used only if the stored point is bit-identical
  to the queried one. `generate_mesh!` drops the points behind the current x slab, so the
  cache holds about one slab.

Gripper case above, the same machine and code otherwise (stages as in the table before; the
warp and slicing rows include their `update_connectivity!`):

| stage | before [s] | calls | after [s] | calls |
|---|---:|---:|---:|---:|
| lattice corners | 1.0 | 0.54 M | 1.0 | 0.54 M |
| A15 fill | 19.5 | 12.70 M | 10.9 | 4.66 M |
| `warp!` (bisection) | 216.1 | 55.90 M | 73.1 | 18.91 M |
| slicing (bisection) | 11.1 | 2.59 M | 2.4 | 0.41 M |
| inverted + isolated + INE | 3.2 | 0 | 3.4 | 0 |
| plane alignment | 9.2 | 0 | 10.4 | 0 |
| **total** | **260.0** | 71.73 M | **101.3** | 24.52 M |

The warp still makes 18.9 M calls, not 11 M: the cache works per oriented edge, and 472,652
oriented crossing edges are bisected (274,644 unordered ones). Bisecting `Xj → Xi` gives
`1 − alpha` only up to the last bits, so one bisection per unordered edge would change the
output. The bitwise check is recorded in [`VERIFICATION.md`](VERIFICATION.md).

## Open

| change | effect (gripper) | output |
|---|---|---|
| `find_surface_nodes` with `face_key` tuples (item 4) | alignment 10.4 s → 1.3 s | different on a few nodes, see below |
| one bisection per unordered edge, `1 − alpha` for the other orientation | warp 18.9 M → ~11 M calls, ~−30 s | last bits of `alpha`; can change which cut point a vertex takes when two are within an ulp |
| Brent/Illinois instead of 40 halvings (smooth field, ~6–8 calls) | warp → ~8 s | same root to ~1e-12, not bit-identical |

Item 4 was tried on 2026-11-03 and not applied, because the output must stay identical.
It changes the result because `warp_mesh_by_planes_sdf!` projects the surface nodes one by
one in the iteration order of the `Set` that `find_surface_nodes` returns, and a projection is
rejected when it would invert a tet next to an already moved node (C-F3). The new key changes
the hash order. Measured with the nodes in ascending order: the alignment of 1 node changes on
the session-5 cantilever and on the session-5 gripper, of 7 and 10 nodes on the gripper test
data (`:bisection`, `:linear`). Most of them end on the plane where the old order left them in
place, 2 the other way. Connectivity and `node_sdf` are unchanged.
