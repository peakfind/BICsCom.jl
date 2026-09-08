module BICsCom

using LinearAlgebra
using Bessels

include("discretization.jl")
export Square, SamplingPoints, get_samplingpoints

include("sphericalDtN.jl")
export CylinderCache, build_cylinder_cache, get_coeff
export BoundaryCache, build_boundary_cache, assemble_dtn

include("boundaryCondition.jl")
export apply_bc
export beta_m, assemble_tbc, apply_tbc!

include("mode.jl")
export ModeField, compute_mode, evaluate_field

include("ext.jl")
export plot_square!, plot_samplingpoints!, plot_min_svals!, plot_field_mode!

end