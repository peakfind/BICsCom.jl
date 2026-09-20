@testset"evaluate_field" begin
    n = 5
    inn = 11.6
    ext = 1.0
    hom = 1.0
    r = 0.3 * 2π
    k = 0.7842
    
    sq = Square([0.0, 0.0], π)
    sp = get_samplingpoints(sq, n)
    
    mf = compute_mode(sp, k, r, inn, ext, hom)

    @test mf.k == k
    @test mf.mode == :te
    @test mf.α == 0.0

    # the field at the sampling points must reproduce A * c by construction
    cydc = build_cylinder_cache(4n, k, r, inn, ext)
    coef = get_coeff(4n, cydc, inn, ext)
    bydc = build_boundary_cache(sp, 4n, k, ext)
    _, A = assemble_dtn(sp, 4n, coef, bydc, k, ext)
    fsp = diag(evaluate_field(mf, sp.cc[1, :], sp.cc[2, :]))
    @test fsp ≈ A * mf.expan_coeffs
    
end

@testset "evaluate_field: Bloch continuation" begin
    n = 10
    inn = 11.56
    r = 0.35 * 2π
    α = 0.0776
    k = 0.4854 # propagating BIC 1 from test_propagating_cylinder.jl

    sq = Square([0.0, 0.0], π)
    sp = get_samplingpoints(sq, n)
    mf = compute_mode(sp, k, r, inn, 1.0, 1.0; α = α)

    @test mf.α == α

    step = 2π / 64
    xs = (-3π + step/2):step:(3π - step/2) # offset by step/2 to avoid fold ties at cell boundaries
    ys = [-0.3, 0.0, 0.5]
    shift = 64
    f =  evaluate_field(mf, sq, xs, ys)
    @test f[shift+1:2*shift, :] ≈ exp(im * α * 2π) .* f[1:shift, :]
    @test f[shift+1:end, :] ≈ exp(im * α * 2π) .* f[1:end-shift, :]
end

@testset"post-processing utilities" begin
    n = 5
    inn = 11.6
    ext = 1.0
    hom = 1.0
    r = 0.3 * 2π
    k = 0.7842
    
    sq = Square([0.0, 0.0], π)
    sp = get_samplingpoints(sq, n)
    
    mf = compute_mode(sp, k, r, inn, ext, hom)

    xs = -π:0.05:π
    ys = -π:0.05:π
    f = evaluate_field(mf, xs, ys)
    normalize_field!(f)
    
    @test maximum(abs, f) ≈ 1.0
end