module Optimization

export optimize_mesh!

using Statistics
using StaticArrays
using LinearAlgebra
using Random

using Implicit2TetMesh.Fundamentals
using Implicit2TetMesh.GenerateMesh

include("OptimizeTetMesh.jl")


end
