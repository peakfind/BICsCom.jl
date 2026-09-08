# functions for compute fields of modes corresponding to BICs

"""
    ModeField

Contain information of modes corresponding to BICs. 

# Fields
- `k::Float64`: the wavenumber
- `inn::Float64`: the dielectric constant inside the cylinders or air-holes
- `ext::Float64`: the dielectric constant outside the cylinders or air-holes
- `r::Float64`: the radius of the cylinders or air-holes
- `orders::Vector{Int64}`: the orders of Bessel functions and Neumann functions
- `coeffs::Matrix{Float64}`: the coefficients in the spherical expansion outside the cylinder, see [`get_coeff`](@ref)
- `expan_coeffs::Vector{ComplexF64}`: the coefficients before special functions
"""
struct ModeField 
    k::Float64 
    inn::Float64
    ext::Float64
    r::Float64
    orders::Vector{Int64}
    coeffs::Matrix{Float64}
    expan_coeffs::Vector{ComplexF64}
end

"""
    compute_mode(sp::SamplingPoints, k, r, inn, ext, homo)

TBW

# Arguments
"""
function compute_mode(sp::SamplingPoints, k, r, inn, ext, homo)
    n = sp.n
    
    odrs = [j - 1 - 2n for j in 1:4n]

    cydc = build_cylinder_cache(4n, k, r, inn, ext)
    coef = get_coeff(4n, cydc, inn, ext)
    bydc = build_boundary_cache(sp, 4n, k, ext)
    Λ, A = assemble_dtn(sp, 4n, coef, bydc, k, ext)
    Δ = apply_bc(Λ, n)
    apply_tbc!(Δ, sp, k; homo = homo)
    
    evals, evecs = eigen(Δ)
    idx = argmin(abs.(evals))
    v = evecs[:, idx]
    
    # Reconstruct the Dirichlet data on 4 edges
    vb = v[1:n]
    vl = v[n+1:2n]
    vt = v[2n+1:end]
    vr = reverse(vl)
    u = [vb; vl; vt; vr]
    
    c = A \ u

    return ModeField(k, inn, ext, r, odrs, coef, c)
end

function evaluate_field(mf::ModeField, xs, ys)
    field = zeros(ComplexF64, length(xs), length(ys))
    ki = mf.k * sqrt(mf.inn)
    ke = mf.k * sqrt(mf.ext)
    trunc = length(mf.orders)
    
    @inbounds for j in eachindex(ys)
        y = ys[j]

        for i in eachindex(xs)
            x = xs[i]
            rd = sqrt(x^2 + y^2)
            θ = atan(y, x)
        
            val = zero(ComplexF64)

            if rd < mf.r 
                for t in 1:trunc 
                    val += mf.expan_coeffs[t] * besselj(mf.orders[t], ki * rd) * exp(im * mf.orders[t] * θ)
                end
            else
                for t in 1:trunc
                    radial = mf.coeffs[1, t] * besselj(mf.orders[t], ke * rd) + mf.coeffs[2, t] * bessely(mf.orders[t], ke * rd)
                    val += mf.expan_coeffs[t] * radial * exp(im * mf.orders[t] * θ)
                end
            end
        
            field[i, j] = val
        end
    end

    return field
end