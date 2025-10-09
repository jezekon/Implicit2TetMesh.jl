"""
    count_negative_determinants(mesh::BlockMesh) -> Int

Iterates through all tetrahedral elements in the mesh, calculates their Jacobian 
determinant, and returns the count of elements with negative determinants.

A negative determinant indicates an inverted (inside-out) element, which is 
typically problematic for finite element analysis.

# Arguments
- `mesh::BlockMesh`: The tetrahedral mesh to analyze

# Returns
- `Int`: Number of elements with negative Jacobian determinants
"""
function count_negative_determinants(mesh::BlockMesh)
    # Initialize counter for elements with negative determinants
    negative_count = 0

    # Iterate through all tetrahedral elements
    for (elem_idx, tet) in enumerate(mesh.IEN)
        # Get the four vertices of the tetrahedron
        vertices = [mesh.X[node_id] for node_id in tet]

        # Calculate edge vectors from the first vertex
        a = vertices[2] - vertices[1]  # Edge vector from vertex 1 to 2
        b = vertices[3] - vertices[1]  # Edge vector from vertex 1 to 3
        c = vertices[4] - vertices[1]  # Edge vector from vertex 1 to 4

        # Calculate Jacobian determinant (proportional to the signed volume)
        # The determinant is computed as the triple scalar product: a·(b×c)
        det_j = dot(a, cross(b, c))

        # Increment counter if determinant is negative
        if det_j < 0
            negative_count += 1
        end
    end

    # Print results
    println("Total elements: $(length(mesh.IEN))")
    println("Elements with negative determinant: $negative_count")

    return negative_count
end
