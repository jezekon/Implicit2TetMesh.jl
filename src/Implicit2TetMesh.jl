module Implicit2TetMesh

# Core data structure and SDF operations:
include("Fundamentals/Fundamentals.jl")
using .Fundamentals

# Generate tetrahedra mesh + warp
include("GenerateMesh/GenerateMesh.jl")
using .GenerateMesh

# Mesh optimization (Laplacian smoothing, quality metrics)
include("Optimization/Optimization.jl")
using .Optimization

# Mesh modification (cutting planes)
include("Modification/Modification.jl")
using .Modification

# Utility functions for mesh analysis
include("Utils/Utils.jl")
using .Utils

# Export types and functions from submodules
export PlaneDefinition

end # module Implicit2TetMesh
