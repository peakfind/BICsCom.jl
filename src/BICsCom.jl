module BICsCom

# Write your package code here.
using LinearAlgebra
using Bessels

include("discretization.jl")
export Square, SamplingPoints, get_samplingpoints

end
