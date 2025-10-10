module Implicit2TetMesh

using JLD2

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

# Add this line to include TetMeshGenerator.jl
include("TetMeshGenerator.jl")

# Exports from module Modification:
export PlaneDefinition, Rectangle, Square, Circle, Ellipse

# Exports from main module:
export MeshGenerationOptions, generate_tetrahedral_mesh

end # module Implicit2TetMesh
