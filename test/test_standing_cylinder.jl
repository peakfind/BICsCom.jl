# Standing waves for dielectric cylinders under TE polarization
# Reference:
#       Z. Hu, Y. Y. Lu, Standing waves on two-dimensional periodic dielectric waveguides, J. Opt. 17 (2015) 065601 (7pp)

@testset "standing wave for dielectric cylinders" begin
    
n = 5
inn = 11.6
ext = 1.0
hom = 1.0

msv = []

# a square with period 2π
sq = Square([0.0, 0.0], π)

# generate sampling points along the boundary of the square
sp = get_samplingpoints(sq, n)

# sweep the frequency
ks = 0.7:1e-5:0.8

for k in ks
    cydc = build_cylinder_cache(4n, k, 0.3*2π, inn, ext)
    coef = get_coeff(4n, cydc, inn, ext)
    bydc = build_boundary_cache(sp, 4n, k, ext)
    Λ, _ = assemble_dtn(sp, 4n, coef, bydc, k, ext)
    Δ = apply_bc(Λ, n)
    apply_tbc!(Δ, sp, k; homo = hom)
    
    s = svdvals(Δ)
    push!(msv, minimum(s))
end

ind = sortperm(msv)
rst = ks[ind[1]]

@test rst ≈ 0.7842 atol=1e-3

end