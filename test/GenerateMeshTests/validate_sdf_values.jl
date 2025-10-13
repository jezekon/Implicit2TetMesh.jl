"""
    validate_node_sdf_values(mesh::BlockMesh, tolerance::Float64=0.05)

Validates that stored SDF values at mesh nodes match the values obtained 
from the implicit SDF function interpolation.

This test function compares `mesh.node_sdf` (stored SDF values) against 
`eval_sdf(mesh, position)` (interpolated from the implicit function) for 
each node in the mesh. Warnings are issued for nodes where the discrepancy 
exceeds the specified tolerance.

# Arguments
- `mesh::BlockMesh`: The tetrahedral mesh containing nodes and SDF values
- `tolerance::Float64=0.05`: Maximum acceptable difference between stored and computed SDF values

# Returns
- `Dict`: Statistics including:
  - `max_error`: Maximum absolute difference found
  - `mean_error`: Mean absolute difference across all nodes
  - `problematic_nodes`: Number of nodes exceeding tolerance
  - `total_nodes`: Total number of nodes checked
  - `error_list`: Vector of (node_index, error) tuples for problematic nodes

"""


# stats = validate_node_sdf_values(mesh, 0.05)

function validate_node_sdf_values(mesh::BlockMesh, tolerance::Float64 = 0.05)
    println("\n" * "="^60)
    println("Validating Node SDF Values Against Implicit Function")
    println("="^60)

    # Initialize error tracking
    total_nodes = length(mesh.X)
    errors = Vector{Float64}(undef, total_nodes)
    problematic_nodes = Vector{Tuple{Int,Float64}}()

    # Check each node
    for i = 1:total_nodes
        # Get stored SDF value
        stored_sdf = mesh.node_sdf[i]

        # Compute SDF value from implicit function interpolation
        computed_sdf = eval_sdf(mesh, mesh.X[i])

        # Calculate absolute error
        error = abs(stored_sdf - computed_sdf)
        errors[i] = error

        # Check if error exceeds tolerance
        if error > tolerance
            push!(problematic_nodes, (i, error))

            # Issue warning for this node
            @warn "Node $i: SDF mismatch exceeds tolerance" stored_sdf computed_sdf error position =
                mesh.X[i]
        end
    end

    # Calculate statistics
    max_error = maximum(errors)
    mean_error = sum(errors) / total_nodes
    std_error = sqrt(sum((errors .- mean_error) .^ 2) / total_nodes)
    num_problematic = length(problematic_nodes)

    # Print summary
    println("\nValidation Summary:")
    println("  Total nodes checked: $total_nodes")
    println("  Tolerance: $tolerance")
    println(
        "  Nodes exceeding tolerance: $num_problematic ($(round(100*num_problematic/total_nodes, digits=2))%)",
    )
    println("  Maximum error: $(round(max_error, digits=6))")
    println("  Mean error: $(round(mean_error, digits=6))")
    println("  Std deviation: $(round(std_error, digits=6))")

    if num_problematic > 0
        println("\n⚠️  $(num_problematic) nodes have SDF mismatches exceeding tolerance!")

        # Show worst 5 cases
        sorted_problems = sort(problematic_nodes, by = x -> x[2], rev = true)
        num_to_show = min(5, length(sorted_problems))
        println("\nWorst $num_to_show cases:")
        for i = 1:num_to_show
            node_idx, err = sorted_problems[i]
            println(
                "  Node $node_idx: error = $(round(err, digits=6)) at position $(mesh.X[node_idx])",
            )
        end
    else
        println("\n✓ All node SDF values are within tolerance!")
    end

    println("="^60 * "\n")

    # Return detailed statistics
    return Dict(
        "max_error" => max_error,
        "mean_error" => mean_error,
        "std_error" => std_error,
        "problematic_nodes" => num_problematic,
        "total_nodes" => total_nodes,
        "error_list" => problematic_nodes,
        "all_errors" => errors,
    )
end

# Info:
# This test should be performed at every step of mesh creation.
