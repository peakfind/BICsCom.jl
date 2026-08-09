module BICsCom

using LinearAlgebra
using Bessels

include("discretization.jl")
export Square, SamplingPoints, get_samplingpoints

include("sphericalDtN.jl")
export CylinderCache, build_cylinder_cache, get_coeff
export BoundaryCache, build_boundary_cache, assemble_dtn

end
