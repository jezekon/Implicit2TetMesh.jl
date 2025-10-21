"""
    remove_isolated_components!(mesh::BlockMesh; keep_largest::Bool=true)

Remove isolated tetrahedra or disconnected regions from the mesh.
Uses graph connectivity analysis to identify and remove elements that don't
share faces with the main body of the mesh.

# Arguments
- `mesh::BlockMesh`: The mesh to process (modified in-place)
- `keep_largest::Bool`: If true, keeps only the largest connected component (default: true)

# Algorithm
1. Build adjacency graph: two tetrahedra are neighbors if they share a face (3 nodes)
2. Find all connected components using BFS
3. Keep only the largest component(s)
4. Update mesh connectivity

# Performance
- Time complexity: O(n_elements + n_faces) - linear in mesh size
- Space complexity: O(n_elements) for visited flags and component tracking

# Returns
- Number of elements removed
"""
function remove_isolated_components!(mesh::BlockMesh; keep_largest::Bool = true)
    n_elements = length(mesh.IEN)
    n_elements == 0 && return 0

    @info "Analyzing mesh connectivity to remove isolated components..."

    # Step 1: Build face-to-elements map for efficient neighbor lookup
    # Key: sorted tuple of 3 node indices (a face)
    # Value: vector of element indices sharing this face
    face_to_elements = Dict{NTuple{3,Int},Vector{Int}}()

    @inbounds for (elem_idx, tet) in enumerate(mesh.IEN)
        # Each tetrahedron has 4 triangular faces
        faces = (
            tuple(sort([tet[1], tet[2], tet[3]])...),  # Face opposite to node 4
            tuple(sort([tet[1], tet[2], tet[4]])...),  # Face opposite to node 3
            tuple(sort([tet[1], tet[3], tet[4]])...),  # Face opposite to node 2
            tuple(sort([tet[2], tet[3], tet[4]])...),   # Face opposite to node 1
        )

        for face in faces
            if haskey(face_to_elements, face)
                push!(face_to_elements[face], elem_idx)
            else
                face_to_elements[face] = [elem_idx]
            end
        end
    end

    # Step 2: Find all connected components using BFS
    visited = falses(n_elements)  # BitArray for fast lookup
    components = Vector{Vector{Int}}()  # Store all components

    for start_idx = 1:n_elements
        visited[start_idx] && continue

        # Start new component with BFS
        component = Int[]
        queue = Int[start_idx]
        visited[start_idx] = true

        while !isempty(queue)
            current = popfirst!(queue)
            push!(component, current)

            # Find neighbors through shared faces
            tet = mesh.IEN[current]
            faces = (
                tuple(sort([tet[1], tet[2], tet[3]])...),
                tuple(sort([tet[1], tet[2], tet[4]])...),
                tuple(sort([tet[1], tet[3], tet[4]])...),
                tuple(sort([tet[2], tet[3], tet[4]])...),
            )

            for face in faces
                # Get all elements sharing this face
                for neighbor_idx in face_to_elements[face]
                    if !visited[neighbor_idx]
                        visited[neighbor_idx] = true
                        push!(queue, neighbor_idx)
                    end
                end
            end
        end

        push!(components, component)
    end

    # Step 3: Identify which components to keep
    n_components = length(components)
    println("  Found $n_components connected component(s)")

    if n_components == 1
        println("  Mesh is fully connected - no isolated components found")
        return 0
    end

    # Print component sizes
    for (i, comp) in enumerate(components)
        println("    Component $i: $(length(comp)) elements")
    end

    # Determine which elements to keep
    if keep_largest
        # Keep only the largest component
        largest_idx = argmax(length.(components))
        elements_to_keep = Set(components[largest_idx])
        println(
            "  Keeping largest component ($largest_idx) with $(length(elements_to_keep)) elements",
        )
    else
        # Keep all components above a threshold (e.g., 1% of largest)
        max_size = maximum(length.(components))
        threshold = max(1, max_size ÷ 100)  # At least 1 element

        elements_to_keep = Set{Int}()
        for comp in components
            if length(comp) >= threshold
                union!(elements_to_keep, comp)
            end
        end
        println(
            "  Keeping $(length(elements_to_keep)) elements from components with ≥$threshold elements",
        )
    end

    # Step 4: Filter mesh elements
    n_original = length(mesh.IEN)
    mesh.IEN = [mesh.IEN[i] for i = 1:n_original if i in elements_to_keep]
    n_removed = n_original - length(mesh.IEN)

    println(
        "  Removed $n_removed isolated elements ($(round(100*n_removed/n_original, digits=2))%)",
    )

    return n_removed
end
