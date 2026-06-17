# ==============================================================================
# Test helpers for Implicit2TetMesh.jl
# ==============================================================================
# Check functions + data-loading conveniences shared by the test suite. This file
# contains NO @testset blocks -- it only defines helpers. It is safe to `include`
# more than once (functions are simply redefined). The check helpers port the
# proven logic from scratch/ (watertight.jl, measure.jl) so the numbers stay
# directly comparable with the recorded regression history.

using Implicit2TetMesh
using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh
using JLD2
using LinearAlgebra
using StaticArrays
using Statistics

# ------------------------------------------------------------------------------
# Paths and data loading (all @__DIR__-based, never pwd-relative)
# ------------------------------------------------------------------------------

# Repository data/ directory, resolved from this file's location (test/helpers/).
const DATA_DIR = normpath(joinpath(@__DIR__, "..", "..", "data"))

# Runtime VTU artifacts land here; the directory is git-ignored and created on demand.
const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "output"))

"""
    ensure_output_dir() -> String

Create `test/output/` if missing and return its absolute path. All tests export
their VTU artifacts here so the tracked test/ tree stays clean.
"""
function ensure_output_dir()
    isdir(OUTPUT_DIR) || mkpath(OUTPUT_DIR)
    return OUTPUT_DIR
end

beam_grid_file() = joinpath(DATA_DIR, "beam", "Z_beam_HEX8_FineGrid_B-1.0_smooth-1.jld2")
beam_sdf_file() = joinpath(DATA_DIR, "beam", "Z_beam_HEX8_FineSDF_B-1.0_smooth-1.jld2")
gripper_grid_file() =
    joinpath(DATA_DIR, "gripper", "Z_robot_gripper_HEX8_FineGrid_B-2.0085_smooth-1.jld2")
gripper_sdf_file() =
    joinpath(DATA_DIR, "gripper", "Z_robot_gripper_HEX8_FineSDF_B-2.0085_smooth-1.jld2")

"""
    load_beam() -> (fine_grid, fine_sdf)

Load the beam grid and SDF arrays from their JLD2 files. Returns the raw arrays
ready for `BlockMesh(fine_sdf, fine_grid)`.
"""
function load_beam()
    gf = beam_grid_file()
    sf = beam_sdf_file()
    @load gf fine_grid
    @load sf fine_sdf
    return (fine_grid, fine_sdf)
end

"""
    load_gripper() -> (fine_grid, fine_sdf)

Load the gripper grid and SDF arrays from their JLD2 files.
"""
function load_gripper()
    gf = gripper_grid_file()
    sf = gripper_sdf_file()
    @load gf fine_grid
    @load sf fine_sdf
    return (fine_grid, fine_sdf)
end

# ------------------------------------------------------------------------------
# Low-level geometry helpers
# ------------------------------------------------------------------------------

"""
    sorted3(a, b, c) -> NTuple{3,Int}

Return the three node indices sorted ascending, as a tuple, without allocating a
temporary array. This is the canonical, orientation-independent key for a
triangular face (same pattern as `face_key` in RemoveIsolatedComponents.jl).
"""
function sorted3(a::Int, b::Int, c::Int)::NTuple{3,Int}
    a > b && ((a, b) = (b, a))
    b > c && ((b, c) = (c, b))
    a > b && ((a, b) = (b, a))
    return (a, b, c)
end

"""
    tet_faces(tet) -> NTuple{4,NTuple{3,Int}}

The four triangular faces of a tetrahedron as sorted vertex tuples.
"""
function tet_faces(tet)
    a, b, c, d = tet[1], tet[2], tet[3], tet[4]
    return (sorted3(a, b, c), sorted3(a, b, d), sorted3(a, c, d), sorted3(b, c, d))
end

# ------------------------------------------------------------------------------
# Invariant checks
# ------------------------------------------------------------------------------

"""
    check_watertight(mesh) -> NamedTuple

Watertightness / manifold check for the trimmed mesh (ported from
scratch/watertight.jl). A tetrahedral mesh has a watertight boundary iff the
boundary faces (incident to exactly ONE tet) form a closed 2-manifold: every
boundary edge is shared by exactly TWO boundary faces.

The single CORRECTNESS criterion is `open_edges == 0` (a crack/hole). The
non-manifold counts are reported for diagnostics only -- they are benign
thin-feature geometry and are NOT invariant across pipeline stages, so the suite
must not assert on them.

Returns `(open_edges, nonmanifold_edges, boundary_faces, nonmanifold_interior)`.
"""
function check_watertight(mesh::BlockMesh)
    # Count how many tets each face belongs to; boundary faces appear exactly once.
    face_count = Dict{NTuple{3,Int},Int}()
    for tet in mesh.IEN
        for f in tet_faces(tet)
            face_count[f] = get(face_count, f, 0) + 1
        end
    end
    boundary_faces = [f for (f, c) in face_count if c == 1]
    nonmanifold_interior = count(c -> c > 2, values(face_count))

    # Every boundary face contributes 3 edges; a closed surface has each edge twice.
    edge_count = Dict{NTuple{2,Int},Int}()
    for f in boundary_faces
        a, b, c = f
        for e in ((min(a, b), max(a, b)), (min(a, c), max(a, c)), (min(b, c), max(b, c)))
            edge_count[e] = get(edge_count, e, 0) + 1
        end
    end
    open_edges = count(c -> c == 1, values(edge_count))        # crack / hole
    nonmanifold_edges = count(c -> c > 2, values(edge_count))

    return (
        open_edges = open_edges,
        nonmanifold_edges = nonmanifold_edges,
        boundary_faces = length(boundary_faces),
        nonmanifold_interior = nonmanifold_interior,
    )
end

"""
    dihedral_stats(mesh) -> NamedTuple

Dihedral-angle statistics + inverted-element count for a mesh, ported unchanged
from scratch/measure.jl `dih()` so the numbers stay comparable with all recorded
history. The dihedral statistics are purely geometric (sign-agnostic).

Returns `(min, max, lt5, lt10, gt140, inverted)`:
  - `min`, `max` : extreme dihedral angles over all tets (degrees)
  - `lt5`, `lt10`: number of dihedral angles below 5 deg / 10 deg
  - `gt140`      : number of dihedral angles above 140 deg
  - `inverted`   : number of tets with non-positive float signed volume
"""
function dihedral_stats(mesh::BlockMesh)
    mn = 180.0
    mx = 0.0
    lt5 = 0
    lt10 = 0
    gt140 = 0
    inverted = 0
    for tet in mesh.IEN
        v = (mesh.X[tet[1]], mesh.X[tet[2]], mesh.X[tet[3]], mesh.X[tet[4]])
        # Float signed-volume sign, identical to scratch/measure.jl for comparability.
        dot(v[2] - v[1], cross(v[3] - v[1], v[4] - v[1])) <= 0 && (inverted += 1)
        # Outward unit normals of the four faces.
        faces = ((1, 2, 3, 4), (1, 2, 4, 3), (1, 3, 4, 2), (2, 3, 4, 1))
        normals = SVector{3,Float64}[]
        ok = true
        for (a, b, c, o) in faces
            cr = cross(v[b] - v[a], v[c] - v[a])
            nlen = norm(cr)
            if nlen < 1e-14
                ok = false
                break
            end
            u = cr / nlen
            ctr = (v[a] + v[b] + v[c]) / 3
            dot(u, v[o] - ctr) > 0 && (u = -u)
            push!(normals, u)
        end
        ok || continue
        for i = 1:3, j = (i+1):4
            ang = acos(-clamp(dot(normals[i], normals[j]), -1.0, 1.0)) * 180 / pi
            mn = min(mn, ang)
            mx = max(mx, ang)
            ang < 5 && (lt5 += 1)
            ang < 10 && (lt10 += 1)
            ang > 140 && (gt140 += 1)
        end
    end
    return (min = mn, max = mx, lt5 = lt5, lt10 = lt10, gt140 = gt140, inverted = inverted)
end

"""
    count_inverted_exact(mesh) -> Int

Number of tets that are NOT strictly positively oriented, decided EXACTLY by the
same predicate the pipeline uses (`is_positively_oriented`). This counts both
inverted (negative volume) and flat (zero volume) tets. On a valid output mesh it
must be 0.

GOTCHA: `ExactPredicates.orient` returns the OPPOSITE sign of the float
determinant; the pipeline helper already encapsulates that, so we route through it
and never call `orient()` raw.
"""
function count_inverted_exact(mesh::BlockMesh)
    n = 0
    for tet in mesh.IEN
        Implicit2TetMesh.GenerateMesh.is_positively_oriented(mesh, tet) || (n += 1)
    end
    return n
end

"""
    min_signed_volume(mesh) -> Float64

Smallest float signed volume over all tets. A positive value means every element
has positive volume (the exact backstop is `count_inverted_exact == 0`). Returns
0.0 for an empty mesh.
"""
function min_signed_volume(mesh::BlockMesh)
    isempty(mesh.IEN) && return 0.0
    vmin = Inf
    for tet in mesh.IEN
        a = mesh.X[tet[2]] - mesh.X[tet[1]]
        b = mesh.X[tet[3]] - mesh.X[tet[1]]
        c = mesh.X[tet[4]] - mesh.X[tet[1]]
        vmin = min(vmin, dot(a, cross(b, c)) / 6.0)
    end
    return vmin
end

"""
    mesh_total_volume(mesh) -> Float64

Sum of the signed volumes of all tets -- the meshed solid volume. On a valid (non-inverted)
mesh every term is positive, so this is the total volume. Used by the Etapa-11 cap-recovery
efficacy check (recovering convex under-cut material must RAISE the total volume). Returns
0.0 for an empty mesh.
"""
function mesh_total_volume(mesh::BlockMesh)
    isempty(mesh.IEN) && return 0.0
    v = 0.0
    for tet in mesh.IEN
        a = mesh.X[tet[2]] - mesh.X[tet[1]]
        b = mesh.X[tet[3]] - mesh.X[tet[1]]
        c = mesh.X[tet[4]] - mesh.X[tet[1]]
        v += dot(a, cross(b, c)) / 6.0
    end
    return v
end

"""
    count_components(mesh) -> Int

Number of connected components, where two tets are connected if they share a face.
Read-only counterpart of `remove_isolated_components!`, reusing the same
face->elements map + BFS pattern from RemoveIsolatedComponents.jl.
"""
function count_components(mesh::BlockMesh)
    n_elements = length(mesh.IEN)
    n_elements == 0 && return 0

    # face -> list of element indices sharing that face
    face_to_elements = Dict{NTuple{3,Int},Vector{Int}}()
    for (elem_idx, tet) in enumerate(mesh.IEN)
        for f in tet_faces(tet)
            push!(get!(face_to_elements, f, Int[]), elem_idx)
        end
    end

    visited = falses(n_elements)
    n_components = 0
    for start_idx = 1:n_elements
        visited[start_idx] && continue
        n_components += 1
        queue = Int[start_idx]
        visited[start_idx] = true
        while !isempty(queue)
            current = popfirst!(queue)
            for f in tet_faces(mesh.IEN[current])
                for neighbor_idx in face_to_elements[f]
                    if !visited[neighbor_idx]
                        visited[neighbor_idx] = true
                        push!(queue, neighbor_idx)
                    end
                end
            end
        end
    end
    return n_components
end

"""
    all_finite_coords(mesh) -> Bool

True iff every node coordinate is finite (no NaN / Inf).
"""
function all_finite_coords(mesh::BlockMesh)
    for p in mesh.X
        all(isfinite, p) || return false
    end
    return true
end

"""
    all_indices_in_bounds(mesh) -> Bool

True iff every element has exactly 4 node indices and all of them lie within
`1:length(mesh.X)`.
"""
function all_indices_in_bounds(mesh::BlockMesh)
    n = length(mesh.X)
    for tet in mesh.IEN
        length(tet) == 4 || return false
        for idx in tet
            (1 <= idx <= n) || return false
        end
    end
    return true
end

"""
    node_sdf_consistency(mesh; pristine_tol=1e-6) -> NamedTuple

Node-SDF consistency invariant for the edge-warp pipeline. There are two node
classes, and ONLY the pristine class admits a tight tolerance:
  - on-surface nodes (`node_sdf == 0`): snapped to the LINEAR cut point by `warp!`
    / slicing, so they intentionally sit OFF the trilinear zero set (`eval_sdf`
    there is non-zero, up to ~0.1 for the beam). These are correct by construction
    and are NOT checked against a tolerance.
  - pristine nodes (`node_sdf != 0`): untouched lattice nodes whose stored value
    must equal the interpolated field `eval_sdf(X)` within `pristine_tol`
    (empirically the error is exactly 0.0 -- they sit on lattice corners).

This is the meaningful, stage-invariant replacement for a blanket `max_error`
bound, which the old Newton-to-isosurface warp satisfied but the edge warp does
not. Returns `(pristine_max_err, n_pristine, n_on_surface, pristine_within_tol)`.
"""
function node_sdf_consistency(mesh::BlockMesh; pristine_tol::Float64 = 1e-6)
    pristine_max_err = 0.0
    n_pristine = 0
    n_on_surface = 0
    for i in eachindex(mesh.X)
        if mesh.node_sdf[i] == 0.0
            n_on_surface += 1
        else
            n_pristine += 1
            err = abs(mesh.node_sdf[i] - eval_sdf(mesh, mesh.X[i]))
            pristine_max_err = max(pristine_max_err, err)
        end
    end
    return (
        pristine_max_err = pristine_max_err,
        n_pristine = n_pristine,
        n_on_surface = n_on_surface,
        pristine_within_tol = pristine_max_err <= pristine_tol,
    )
end

"""
    boundary_max_abs_sdf(mesh) -> Float64

Largest `|eval_sdf|` over all BOUNDARY vertices (vertices of faces incident to exactly one
tet) -- i.e. how far the mesh boundary strays from the implicit surface, measured in field
units. This is the surface-fidelity invariant of the `cut_points = :bisection` mode: every
boundary vertex then sits on the trilinear zero set, so the value is ~0 (vs ~0.13 on the
beam with `:linear` on the non-distance smoothed input). Returns 0.0 for an empty mesh.
"""
function boundary_max_abs_sdf(mesh::BlockMesh)
    face_count = Dict{NTuple{3,Int},Int}()
    for tet in mesh.IEN
        for f in tet_faces(tet)
            face_count[f] = get(face_count, f, 0) + 1
        end
    end
    boundary_verts = Set{Int}()
    for (f, c) in face_count
        if c == 1
            push!(boundary_verts, f[1])
            push!(boundary_verts, f[2])
            push!(boundary_verts, f[3])
        end
    end
    m = 0.0
    for v in boundary_verts
        m = max(m, abs(eval_sdf(mesh, mesh.X[v])))
    end
    return m
end

"""
    boundary_vertices(mesh) -> Set{Int}

Set of vertices that lie on the mesh boundary (vertices of faces incident to exactly one
tetrahedron). Shared by the relaxation invariant checks below.
"""
function boundary_vertices(mesh::BlockMesh)
    face_count = Dict{NTuple{3,Int},Int}()
    for tet in mesh.IEN
        for f in tet_faces(tet)
            face_count[f] = get(face_count, f, 0) + 1
        end
    end
    verts = Set{Int}()
    for (f, c) in face_count
        if c == 1
            push!(verts, f[1]); push!(verts, f[2]); push!(verts, f[3])
        end
    end
    return verts
end

"""
    max_interior_sdf(mesh) -> Float64

Largest `eval_sdf` over the INTERIOR (non-boundary) vertices. The relaxation gate keeps
every moved interior vertex strictly inside, and lattice interior vertices start inside,
so on a valid relaxed mesh this is `< 0`. Returns `-Inf` if there are no interior
vertices.
"""
function max_interior_sdf(mesh::BlockMesh)
    bverts = boundary_vertices(mesh)
    m = -Inf
    for v in eachindex(mesh.X)
        v in bverts && continue
        m = max(m, eval_sdf(mesh, mesh.X[v]))
    end
    return m
end

"""
    surface_edge_cv(mesh) -> Float64

Coefficient of variation (std / mean) of the lengths of all edges incident to a boundary
vertex. The `:uniform` relaxation mode equalizes element sizes near the boundary, so this
must DECREASE across the pass. Returns 0.0 if there are no such edges.
"""
function surface_edge_cv(mesh::BlockMesh)
    bverts = boundary_vertices(mesh)
    seen = Set{NTuple{2,Int}}()
    lengths = Float64[]
    for tet in mesh.IEN
        a, b, c, d = tet[1], tet[2], tet[3], tet[4]
        for (u, w) in ((a, b), (a, c), (a, d), (b, c), (b, d), (c, d))
            (u in bverts || w in bverts) || continue
            ek = (min(u, w), max(u, w))
            ek in seen && continue
            push!(seen, ek)
            push!(lengths, norm(mesh.X[u] - mesh.X[w]))
        end
    end
    isempty(lengths) && return 0.0
    return std(lengths) / mean(lengths)
end

"""
    validate_node_sdf_values(mesh, tolerance=0.005) -> Dict

Silent port of the former test/GenerateMeshTests/validate_sdf_values.jl. Compares
the stored `mesh.node_sdf` against `eval_sdf(mesh, position)` (the implicit field)
at every node and returns the statistics Dict only -- no printing, no warnings.

Keys: `"max_error"`, `"mean_error"`, `"problematic_nodes"` (count over tolerance),
`"total_nodes"`, `"error_list"` (vector of `(node_index, error)` over tolerance).
"""
function validate_node_sdf_values(mesh::BlockMesh, tolerance::Float64 = 0.005)
    total_nodes = length(mesh.X)
    errors = Vector{Float64}(undef, total_nodes)
    problematic_nodes = Tuple{Int,Float64}[]
    for i = 1:total_nodes
        stored_sdf = mesh.node_sdf[i]
        computed_sdf = eval_sdf(mesh, mesh.X[i])
        err = abs(stored_sdf - computed_sdf)
        errors[i] = err
        err > tolerance && push!(problematic_nodes, (i, err))
    end
    max_error = total_nodes == 0 ? 0.0 : maximum(errors)
    mean_error = total_nodes == 0 ? 0.0 : sum(errors) / total_nodes
    return Dict(
        "max_error" => max_error,
        "mean_error" => mean_error,
        "problematic_nodes" => length(problematic_nodes),
        "total_nodes" => total_nodes,
        "error_list" => problematic_nodes,
    )
end
