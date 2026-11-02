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

## Bitwise identity of the merge (2026-11-02)

Final mesh (`X`, `IEN`, compared with `isequal`) of beam / gripper (`:linear`, `:bisection`,
with planes), the experimental passes, and the thesis stage-3 production runs (session-5 clipped
RBF fields of the cantilever and the gripper):

- `dev` with `3b099e0` reverted = tag `thesis-recompute-s5` in all 8 cases;
- `dev` after the clean-up commits = the merge commit `6a77ab7` in all 14 cases;
- the merge with C-F1 and C-F3 reverted = `origin/dev` before the merge (`41a4378`) in all 13 cases.

## Known issue found by the comparison

`BlockMesh(fine_sdf, fine_grid)` uses the input grid itself as the lattice. When the grid is
stored in Float32 with a spacing that is not exact in binary (gripper test data: h = 2.0085,
nodes up to 4e-5 off a uniform grid), a tile vertex shared by two cells gets slightly different
coordinates from each, and the node hash (`quantize` with tolerance 1e-8 h) does not merge them:
the fill has 480 682 nodes instead of 401 876. The copies are warped independently before
`merge_duplicate_nodes!` runs, and the final gripper mesh keeps 965 non-manifold boundary edges
and 43 node pairs closer than 1e-3 h; the same field on an exactly uniform lattice
(`BlockMesh(source; dx)`) gives 0 and 0. quartet keys its vertices by integer lattice
coordinates and cannot have this. Not affected: the beam data (h = 1, exact), the
`BlockMesh(source; dx)` path (lattice built in Float64), and so the thesis-recompute meshes.
Fix: key the lattice vertices by their integer tile coordinates, as quartet does.
