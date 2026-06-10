# ----------------------------
# Quantization function – unchanged
# ----------------------------
function quantize(p::SVector{3,Float64}, tol::Float64)
    return (round(p[1] / tol) * tol, round(p[2] / tol) * tol, round(p[3] / tol) * tol)
end

# ----------------------------
# Function for discretizing a cell using A15 scheme
# ----------------------------
function process_cell_A15!(mesh::BlockMesh, i::Int, j::Int, k::Int)
    tol = mesh.grid_tol

    # First check the SDF values of the current cell
    current_sdf_values = get_cell_sdf_values(mesh, i, j, k)

    # If any value in the current cell is negative or close to zero (inside or on the
    # surface), we definitely process the cell.
    if any(x -> x <= tol, current_sdf_values)
        # Continue with regular processing
    else
        # All values in the current cell are positive (fully outside) - check neighboring cells

        # Define offsets for neighboring cells (direct neighbors in all directions)
        neighbor_offsets = [
            (1, 0, 0),
            (-1, 0, 0),   # neighbors in x direction
            (0, 1, 0),
            (0, -1, 0),   # neighbors in y direction
            (0, 0, 1),
            (0, 0, -1),    # neighbors in z direction
        ] #TODO: select only relevant one

        # Initialize a flag indicating that all neighboring cells are fully outside
        # (all SDF values positive under the phi < 0 = inside convention).
        all_neighbors_outside = true

        # Check SDF values of neighboring cells
        for (di, dj, dk) in neighbor_offsets
            # Calculate indices of neighboring cell
            ni, nj, nk = i + di, j + dj, k + dk

            # Check grid boundaries
            if 1 <= ni < mesh.nx && 1 <= nj < mesh.ny && 1 <= nk < mesh.nz
                # Get SDF values for neighboring cell
                neighbor_sdf = get_cell_sdf_values(mesh, ni, nj, nk)

                # If any value in the neighboring cell is negative or zero (inside or on
                # the surface), clear the flag and end the check
                if any(x -> x <= tol, neighbor_sdf)
                    all_neighbors_outside = false
                    break
                end
            end
        end

        # If the current cell and all neighboring cells are fully outside,
        # we can safely skip the cell
        if all_neighbors_outside
            return
        end
    end

    # Retrieve min and max corners of the cell
    v000 = mesh.grid[i, j, k]
    v111 = mesh.grid[i+1, j+1, k+1]
    vmins = v000
    vmaxs = v111

    # Precompute differences for coordinate interpolation
    Δ = vmaxs .- vmins
    local_mapping = Dict{Int,Int}()

    @inbounds for li = 1:length(tile_ref)
        # tile_ref is assumed to be defined in A15_scheme.jl and normalized (in [0,1]*4 originally)
        local_coord = SVector{3,Float64}(tile_ref[li] ./ 4.0)  # normalized coordinates in [0,1]
        # Compute physical point using linear interpolation
        p = vmins .+ Δ .* local_coord
        p = SVector{3,Float64}(p)  # ensure static vector type
        p_key = quantize(p, tol)
        if haskey(mesh.node_hash, p_key)
            local_mapping[li] = mesh.node_hash[p_key]
        else
            push!(mesh.X, p)

            sdf_aprox = eval_sdf(mesh, p)
            push!(mesh.node_sdf, sdf_aprox)
            local_index = length(mesh.X)
            local_mapping[li] = local_index
            mesh.node_hash[p_key] = local_index
        end
    end

    # Process tetrahedral connectivity from A15 scheme
    @inbounds for tet in tetra_connectivity
        global_tet = [local_mapping[li] for li in tet]

        # Get the coordinates of the tetrahedron vertices
        tet_coords = [mesh.X[idx] for idx in global_tet]

        # Directly evaluate SDF at each vertex position for maximum accuracy
        # This is more accurate than using pre-computed values
        tet_sdf = [eval_sdf(mesh, coord) for coord in tet_coords]

        # Update the stored SDF values with these more accurate evaluations
        for (i, idx) in enumerate(global_tet)
            mesh.node_sdf[idx] = tet_sdf[i]
        end

        # Include tetrahedron only if at least one vertex is inside or on the boundary
        if any(x -> x <= 0, tet_sdf)
            push!(mesh.IEN, global_tet)
        end
    end
end

# ----------------------------
# Merge duplicate nodes after mesh generation (unchanged logic, only type annotations updated)
# ----------------------------
function merge_duplicate_nodes!(mesh::BlockMesh)
    @info "Merging duplicate nodes"
    println("  Before merging duplicates: $(length(mesh.X)) nodes")
    tol = mesh.grid_tol*10^4
    new_nodes = Vector{SVector{3,Float64}}()
    new_node_sdf = Vector{Float64}()
    node_map = Dict{Int,Int}()
    global_hash = Dict{NTuple{3,Float64},Int}()
    @inbounds for i = 1:length(mesh.X)
        p = mesh.X[i]
        p_key = quantize(p, tol)
        if haskey(global_hash, p_key)
            node_map[i] = global_hash[p_key]
        else
            push!(new_nodes, p)
            push!(new_node_sdf, mesh.node_sdf[i])
            new_index = length(new_nodes)
            node_map[i] = new_index
            global_hash[p_key] = new_index
        end
    end
    @inbounds for i = 1:length(mesh.IEN)
        mesh.IEN[i] = [node_map[old] for old in mesh.IEN[i]]
    end
    mesh.X = new_nodes
    mesh.node_sdf = new_node_sdf
    mesh.node_map = node_map
    println("  After merging duplicates:  $(length(mesh.X)) nodes")
end

# ----------------------------
# Cleanup unused nodes and reindex connectivity (unchanged logic)
# ----------------------------
function cleanup_unused_nodes!(mesh::BlockMesh)
    @info "Cleaning unused nodes"
    println("  Number of nodes before cleanup: $(length(mesh.X))")
    used_nodes = Set{Int64}()
    @inbounds for element in mesh.IEN
        union!(used_nodes, element)
    end
    new_node_map = Dict{Int64,Int64}()
    new_coords = Vector{SVector{3,Float64}}()
    new_node_sdf = Vector{Float64}()
    sorted_used = sort(collect(used_nodes))
    for (new_id, old_id) in enumerate(sorted_used)
        new_node_map[old_id] = new_id
        push!(new_coords, mesh.X[old_id])
        push!(new_node_sdf, mesh.node_sdf[old_id])
    end
    @inbounds for i = 1:length(mesh.IEN)
        mesh.IEN[i] = [new_node_map[old_id] for old_id in mesh.IEN[i]]
    end
    mesh.X = new_coords
    mesh.node_sdf = new_node_sdf
    mesh.node_map = new_node_map
    println("  Number of nodes after cleanup:  $(length(mesh.X))")
end

# ----------------------------
# Create inverse connectivity (unchanged logic)
# ----------------------------
function create_INE!(mesh::BlockMesh)
    mesh.INE = [Vector{Int64}() for _ = 1:length(mesh.X)]
    @inbounds for (elem_id, element) in enumerate(mesh.IEN)
        for node_id in element
            push!(mesh.INE[node_id], elem_id)
        end
    end
    return mesh
end

# ----------------------------
# Modified mesh generation function with scheme selection
# ----------------------------
function generate_mesh!(mesh::BlockMesh, scheme::String)
    empty!(mesh.X)
    empty!(mesh.IEN)
    empty!(mesh.node_sdf)
    empty!(mesh.node_map)
    empty!(mesh.cell_center_map)
    empty!(mesh.node_hash)

    @info "Generating mesh with $(scheme) scheme..."

    @inbounds for i = 1:(mesh.nx-1)
        for j = 1:(mesh.ny-1)
            for k = 1:(mesh.nz-1)
                if scheme == "A15"
                    process_cell_A15!(mesh, i, j, k)
                else
                    error("Unknown scheme: $scheme. Only 'A15' is supported.")
                end
            end
        end
    end

    cleanup_unused_nodes!(mesh)
    create_INE!(mesh)
end


# ------------------------------------------------------------------
# Edge-based warp to linear cut points (port of quartet warp_vertices)
# ------------------------------------------------------------------
"""
    warp!(mesh, scheme; threshold = 0.3)

Snap lattice vertices onto the linear cut points of incident sign-crossing edges.

For every tetrahedron edge `(i, j)` whose endpoints have strictly opposite SDF signs the
isosurface crosses the edge at the linear cut point

    X_cut = X_i + alpha * (X_j - X_i),   alpha = phi_i / (phi_i - phi_j)  in (0, 1).

If that cut point is close to an endpoint (`alpha < threshold` for `i`, or
`alpha > 1 - threshold` for `j`), the endpoint is moved onto the cut point and its stored
SDF is set to zero, i.e. the vertex is placed exactly on the surface. Each vertex is
warped to its CLOSEST qualifying cut point (smallest `alpha * |edge|^2`, quartet's metric).
Displacements are computed from the original positions and SDF values and applied all at
once, so the result does not depend on vertex order.

The vertex is moved to the LINEAR cut point only (no Newton projection to the true
isosurface); keeping the move linear is what leaves the surrounding tetrahedra well shaped.
This is a direct port of quartet's `warp_vertices` (make_tet_mesh.cpp:152-193); see
Labelle 2007 §3.2.

`threshold` is quartet's warp coefficient (must lie in `[0, 0.5]`) and is exposed so later
stages can tune it. `scheme` is kept only for call-site compatibility: the warp operates on
`mesh.IEN` / `mesh.node_sdf` and is identical for every scheme.
"""
function warp!(mesh::BlockMesh, scheme::String; threshold::Float64 = 0.3)
    @info "Warping vertices to edge cut points (threshold = $threshold)..."
    @assert 0.0 <= threshold <= 0.5 "warp threshold must lie in [0, 0.5]"

    n = length(mesh.X)
    # Best (closest) qualifying cut point found so far, per vertex.
    best_metric = fill(Inf, n)                              # quartet's warp[] : alpha * |edge|^2
    displacement = [zero(SVector{3,Float64}) for _ = 1:n]   # quartet's d[]    : vector onto the cut point
    to_warp = falses(n)                                     # quartet's warp_nbr >= 0

    # Visit every edge of every tetrahedron. Shared edges are visited several times,
    # but keeping the per-vertex minimum makes the repetition harmless.
    for tet in mesh.IEN
        for u = 1:3
            i = tet[u]
            for v = (u+1):4
                j = tet[v]
                phi_i = mesh.node_sdf[i]
                phi_j = mesh.node_sdf[j]

                # Only strictly sign-crossing edges have an interior cut point.
                if (phi_i < 0 && phi_j > 0) || (phi_i > 0 && phi_j < 0)
                    # Fraction along the edge from i to the zero crossing, in (0, 1).
                    alpha = phi_i / (phi_i - phi_j)
                    edge2 = sum(abs2, mesh.X[j] - mesh.X[i])

                    if alpha < threshold
                        # Cut point is close to i -> warp i toward j.
                        metric = alpha * edge2
                        if metric < best_metric[i]
                            best_metric[i] = metric
                            displacement[i] = alpha * (mesh.X[j] - mesh.X[i])
                            to_warp[i] = true
                        end
                    elseif alpha > 1 - threshold
                        # Cut point is close to j -> warp j toward i.
                        metric = (1 - alpha) * edge2
                        if metric < best_metric[j]
                            best_metric[j] = metric
                            displacement[j] = (1 - alpha) * (mesh.X[i] - mesh.X[j])
                            to_warp[j] = true
                        end
                    end
                end
            end
        end
    end

    # Apply all warps at once: move the vertex onto the surface and mark it as on it.
    warped_count = 0
    for v = 1:n
        if to_warp[v]
            mesh.X[v] = mesh.X[v] + displacement[v]
            mesh.node_sdf[v] = 0.0
            warped_count += 1
        end
    end
    println("  Warped $warped_count vertices onto the isosurface")
end

# ---------------------------------------------------
# Function: Update mesh topology (mesh.X, mesh.IEN, mesh.INE)
# ---------------------------------------------------
function update_connectivity!(mesh::BlockMesh)
    cleanup_unused_nodes!(mesh)        # Recalculates mesh.X, mesh.node_sdf and reindexes mesh.IEN and mesh.node_map
    merge_duplicate_nodes!(mesh)       # Merges duplicate nodes and adjusts connectivity in mesh.IEN
    create_INE!(mesh)                  # Creates inverse connectivity (mesh.INE)
end
