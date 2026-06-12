module Fundamentals

export BlockMesh, get_cell_sdf_values, eval_sdf, calculate_volume_from_sdf
export SDFSource, StructuredSDF, UnstructuredSDF, bbox
export build_sdf_source, simp_to_sdf_source, levelset_to_sdf_source

using StaticArrays
using LinearAlgebra
using FastGaussQuadrature
using Base.Threads

# Pluggable SDF sources (Etapa 8). Included before BlockMesh.jl, which holds an
# SDFSource and builds a StructuredSDF for the structured path. Order matters:
# Adapters.jl uses hex_volume from UnstructuredSDF.jl, which uses the shape
# functions; SDFSource.jl declares the abstract type both concrete sources extend.
include("SDFSources/SDFSource.jl")
include("SDFSources/StructuredSDF.jl")
include("SDFSources/UnstructuredSDF.jl")
include("SDFSources/Adapters.jl")

include("BlockMesh.jl")
include("SDFOperations.jl")
include("CalcVolumeFromSDF.jl")

end
