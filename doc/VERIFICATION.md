# Verification against quartet

quartet (Doran, Chang & Bridson 2013, <https://github.com/crawforddoran/quartet>) is the
reference implementation of isosurface stuffing on the A15 lattice. The package's lattice
fill, warp and trimming are ports of its `cut_lattice`, `warp_vertices` and `trim_spikes`;
the quadruple-zero pass extends its `remove_exterior_tets` with a shape filter and a pocket
rule. This file records how the two are compared and the last result.

## Method (`test/quartet_crossval/`)

- **Same input.** quartet reads an SDF grid and places one A15 tile on 4×4×4 grid cells, with
  its vertices on grid nodes; the package places one tile on every cell of its lattice, with
  vertices at multiples of h/4. On the input grid refined 4× the two lattices coincide. The
  package's trilinear field is sampled on that refined grid and rounded to Float32 (quartet's
  precision); quartet gets these values, the package meshes a `StructuredSDF` of the same
  values on a lattice of spacing h.
- **Unmodified quartet.** `quartet_driver.cpp` includes `make_tet_mesh.cpp` as a source file and
  calls its static functions in the order of `make_tet_mesh` (`optimize = false`, no features),
  dumping the mesh after each step in full Float32 precision.
- **Step by step.** `xval.jl` runs the package steps (`generate_mesh!`, `warp!`,
  `update_connectivity!`, `slice_ambiguous_tetrahedra!`) and compares vertices by position and
  tetrahedra as vertex sets after each step. Differences in the quadruple-zero tets are
  classified by the two documented package rules (dihedral filter [10°, 140°], pocket rule).
- Not part of `Pkg.test()` (needs a C++ build); how to build and run is in the header of
  `xval.jl`.

## Result (2026-11-02, the code after the merge of the thesis-recompute fixes into `dev`)

| step | beam (h = 1) | gripper (h = 2.0085) |
|---|---|---|
| lattice (`cut_lattice` vs `generate_mesh!`) | identical: 30 510 vertices, 135 401 tets | identical: 401 876 vertices, 2 056 474 tets |
| warp | same 7 614 warped vertices; positions within 2.1e-6 h (Float32) | same 48 264 warped vertices; positions within 1.1e-5 h, except 2 vertices (below) |
| trim + quadruple-zero (vs `trim_spikes` + `remove_exterior_tets`) | identical: 92 679 tets, of them 735 quadruple-zero | identical: 1 982 554 tets, of them 86 quadruple-zero, except 27 tets around the 2 vertices |

- The 2 gripper vertices have two candidate cut points whose warp metrics `α·|e|²` differ by
  3e-6 (relative). The package (Float64) takes the smaller; quartet's Float32 positions round the
  order the other way. A precision effect of quartet, not a difference of the rule.
- On these inputs the shape filter and the pocket rule of the quadruple-zero pass change nothing
  against quartet's centroid test.
- After the slice the package also removes face-disconnected pieces
  (`remove_isolated_components!`): beam 4 tets, gripper 0.
- The same comparison with the quadruple-zero rule as it was before `3b099e0` (an extra
  "no face adjoins the solid → discard" rule) gives 147 well-shaped quadruple-zero tets on the
  beam that quartet keeps and the package discarded; `3b099e0` removed that rule.
- Rerun after the integer lattice keys (below): the same result in every row.

## Bitwise identity of the merge (2026-11-02)

Final mesh (`X`, `IEN`, compared with `isequal`) of beam / gripper (`:linear`, `:bisection`,
with planes), the experimental passes, and the thesis stage-3 production runs (session-5 clipped
RBF fields of the cantilever and the gripper):

- `dev` with `3b099e0` reverted = tag `thesis-recompute-s5` in all 8 cases;
- `dev` after the clean-up commits = the merge commit `6a77ab7` in all 14 cases;
- the merge with C-F1 and C-F3 reverted = `origin/dev` before the merge (`41a4378`) in all 13 cases.

## Bitwise identity of the evaluation caches (2026-11-03)

`CutAlphaCache` (warp and slicing) and `TileSDFCache` (A15 fill), see
[`PERFORMANCE_TODO.md`](PERFORMANCE_TODO.md). Final mesh (`X`, `node_sdf`, `IEN` with
`isequal`) against `2e51c5a`: identical in all 6 cases, i.e. the thesis stage-3 pipeline on the
session-5 clipped fields of the cantilever and the gripper (`:bisection`, with the alignment
planes), and `generate_tetrahedral_mesh` with planes on the beam and on the gripper test data
(Float32 grid, not exactly uniform), each with `:linear` and `:bisection`.
`test/runtests.jl` checks that a shared cache gives the same mesh as a separate cache per stage.

## Issue found by the comparison (fixed)

`BlockMesh(fine_sdf, fine_grid)` uses the input grid itself as the lattice. When the grid is
stored in Float32 with a spacing that is not exact in binary (gripper test data: h = 2.0085,
nodes up to 4e-5 off a uniform grid; the bunny data too), a tile vertex shared by two cells got
slightly different coordinates from each, and the node hash (position rounded to 1e-8 h) did
not merge them: the gripper fill had 480 682 nodes instead of 401 876. The copies were warped
independently before `merge_duplicate_nodes!` runs; where they ended in different places, the
neighbouring tets no longer shared the vertex and the mesh was cracked inside. The cracks are
closed surfaces, so the watertightness check (open edges) did not see them.

Fix: `process_cell_A15!` keys the tile vertices by their integer lattice coordinate in quarter
cells, as quartet does (`Vec3i`). Regression test: `test/test_nonuniform_grid.jl`.

Gripper test data through the file path, default options:

| | before | after | same field, uniform Float64 lattice |
|---|---|---|---|
| nodes / tets | 414 072 / 1 984 321 | 408 614 / 1 982 554 | 408 614 / 1 982 554 |
| boundary faces | 273 246 | 238 346 | 238 346 |
| of them with the centroid deeper than 0.05 h inside | 34 510 | 1 127 | 1 127 |
| non-manifold boundary edges | 965 | 0 | 0 |
| node pairs closer than 1e-3 h | 43 | 0 | 0 |

The 1 127 deep faces left are the sagitta of curved surface faces, not cracks; the tet count is
quartet's (table above).

Unchanged bitwise (`X`, `node_sdf`, `IEN` with `isequal`): the fill of the beam data (h = 1,
exact in binary) and of every `BlockMesh(source; dx)` lattice (built in Float64): the thesis
session-5 cantilever and gripper fields and the bunny benchmark at dx = 4, 6, 8; the final mesh
of the beam (`:linear`, `:bisection`) and of the session-5 cantilever production run. The code
after the fill did not change, so an identical fill gives an identical mesh. Not affected
therefore: the thesis-recompute meshes and the bunny benchmark table.
