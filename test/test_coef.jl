@testset "Computing the coeffcients" begin
    # assume that the medium is homogeneous
    k = 0.5
    n = 5
    inn = 1.0
    ext = 1.0

    # generate the coeffcients
    sq = Square([0.0, 0.0], π) # the period is 2π
    sp = get_samplingpoints(sq, n)
    cydc = build_cylinder_cache(4n, k, 0.3*2π, inn, ext) # the radius is arbitrary
    coef = get_coeff(4n, cydc, inn, ext)
    
    # exact coefficients
    exact_coef = Matrix{Float64}(undef, 2, 4n)
    exact_coef[1, :] .= 1.0
    exact_coef[2, :] .= 0.0
    
    @test coef ≈ exact_coef
end