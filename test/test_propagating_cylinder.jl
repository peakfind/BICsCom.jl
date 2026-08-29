# Propagating waves for dielectric cylinders under TE polarization
# Reference:
#     L. Yuan, Y. Y. Lu, Propagating Bloch modes above the lightline on a periodic array of cylinders, J. Phys. B: At. Mol. Opt. Phys. 50 (2017) 05LT01 (5pp) 
#
# ϵ₁ = 11.56, r = 0.35 * 2π
# propagating BIC 1: x-even, βL/2π = 0.0776, ωL/2πc = 0.4854
# propagating BIC 2: x-odd BIC: βL/2π = 0.2483, ωL/2πc = 0.6702

@testset"propagating wave for dielectric cylinders" begin
    n = 10
    inn = 11.56
    ext = 1.0
    hom = 1.0
    
    @testset"propagating BIC 1" begin
        α = 0.0776
        msv = []
        sq = Square([0.0, 0.0], π)
        sp = get_samplingpoints(sq, n)
        ks = 0.4:1e-5:0.5

        for k in ks
            cydc = build_cylinder_cache(4n, k, 0.35*2π, inn, ext)
            coef = get_coeff(4n, cydc, inn, ext)
            bydc = build_boundary_cache(sp, 4n, k, ext)
            Λ, _ = assemble_dtn(sp, 4n, coef, bydc, k, ext)
            Δ = apply_bc(Λ, n; α = α)
            apply_tbc!(Δ, sp, k; α = α, homo = hom)
            
            s = svdvals(Δ)
            push!(msv, minimum(s))
        end
        
        ind = sortperm(msv)
        rst = ks[ind[1]]
        @test rst ≈ 0.4854 atol=1e-3
    end
    
    @testset"propagating BIC 2" begin
        α = 0.2483
        msv = []
        sq = Square([0.0, 0.0], π)
        sp = get_samplingpoints(sq, n)
        ks = 0.6:1e-5:0.7

        for k in ks
            cydc = build_cylinder_cache(4n, k, 0.35*2π, inn, ext)
            coef = get_coeff(4n, cydc, inn, ext)
            bydc = build_boundary_cache(sp, 4n, k, ext)
            Λ, _ = assemble_dtn(sp, 4n, coef, bydc, k, ext)
            Δ = apply_bc(Λ, n; α = α)
            apply_tbc!(Δ, sp, k; α = α, homo = hom)
            
            s = svdvals(Δ)
            push!(msv, minimum(s))
        end

        ind = sortperm(msv)
        rst = ks[ind[1]]
        @test rst ≈ 0.6702 atol=1e-3
    end
end