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


end # module Implicit2TetMesh
