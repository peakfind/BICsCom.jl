# Some function to impose the constraints from:
# - periodic bc (standing)
# - quasi-periodic bc (propagating)
# - TBC

"""
    apply_bc(Λ, n::Int64; α::Float64 = 0.0, period::Float64 = 2π)

Impose the periodic boundary condition or quasi-period boundary condition for standing Bloch modes 
or propagating Bloch modes, respectively.

# Arguments
- `Λ`: the DtN matrix obtained by the spherical expansion
- `n`: the number of sampling points on each edge
- `α`: the Bloch wavenumber
- `period`: the period of the periodic slab
"""
function apply_bc(Λ, n::Int64; α::Float64 = 0.0, period::Float64 = 2π)
    
    Δ = Λ[1:3n, 1:3n]

    if α == 0.0 # for standing Bloch modes
        @views Δ[1:n, n+1:2n] .+= reverse(Λ[1:n, 3n+1:4n], dims=2)                                         
        @views Δ[2n+1:3n, n+1:2n] .+= reverse(Λ[2n+1:3n, 3n+1:4n], dims=2)                                     
                                                                                                               
        @views Δ[n+1:2n, 1:n] .+= reverse(Λ[3n+1:4n, 1:n], dims=1)                                         
        @views Δ[n+1:2n, 2n+1:3n] .+= reverse(Λ[3n+1:4n, 2n+1:3n], dims=1)                                     
                                                                                                               
        @views Δ[n+1:2n, n+1:2n] .+= reverse(Λ[n+1:2n, 3n+1:4n], dims=2)                                      
        @views Δ[n+1:2n, n+1:2n] .+= reverse(Λ[3n+1:4n, n+1:2n], dims=1)                                      
        @views Δ[n+1:2n, n+1:2n] .+= reverse(Λ[3n+1:4n, 3n+1:4n], dims=:)
    else # for propagating Bloch modes
        e⁺ = exp(im * α * period)                                                                                 
        e⁻ = exp(-im * α * period)                                                                                 
                                                                                                               
        @views Δ[1:n, n+1:2n] .+= e⁺ .* reverse(Λ[1:n, 3n+1:4n], dims=2)                                   
        @views Δ[2n+1:3n, n+1:2n] .+= e⁺ .* reverse(Λ[2n+1:3n, 3n+1:4n], dims=2)                               
                                                                                                               
        @views Δ[n+1:2n, 1:n] .+= e⁻ .* reverse(Λ[3n+1:4n, 1:n], dims=1)                                   
        @views Δ[n+1:2n, 2n+1:3n] .+= e⁻ .* reverse(Λ[3n+1:4n, 2n+1:3n], dims=1)                               
                                                                                                               
        @views Δ[n+1:2n, n+1:2n] .+= e⁺ .* reverse(Λ[n+1:2n, 3n+1:4n], dims=2)                                
        @views Δ[n+1:2n, n+1:2n] .+= e⁻ .* reverse(Λ[3n+1:4n, n+1:2n], dims=1)                                
        @views Δ[n+1:2n, n+1:2n] .+= reverse(Λ[3n+1:4n, 3n+1:4n], dims=:)  
    end
    
    return Δ
end

"""
    beta_m(m::Int64, k; α = 0.0, homo = 1.0, period = 2π)

Compute ``\\beta_{m}`` given by 
```math
\\beta_m = \\sqrt{k^2 - (m + \\frac{2\\pi}{\\text{period}})^2}.
```

# Arguments
- `m`: the order of the Fourier mode
- `k`: the wavenumber
- `α`: the Bloch wavenumber
- `homo`: the dielectric constant of the homogeneous medium outside the periodic slab
- `period`: the period of the periodic slab
"""
function beta_m(m::Int64, k; α = 0.0, homo = 1.0, period = 2π)
    α² = (α + 2π * m / period)^2
    k² = k * k * homo
    
    if k² > α²
        β = complex(sqrt(k² - α²))
    else
        β = im * sqrt(α² - k²)
    end
    
    return β
end


"""
    assemble_tbc(x, beta, fmodes; α = 0.0, period = 2π)

Construct the Transparent Boundary Condition matrix for top and bottom boundaries.

# Arguments
- `x`: x coordinates of sampling points on the top or the bottom boundary
- `beta`: get from [`beta_m`](@ref).
- `fmodes`: the orders in the truncated Fourier serie
- `α`: the Bloch wavenumber
- `period`: the period of the periodic slab
"""
function assemble_tbc(x, beta, fmodes; α = 0.0, period = 2π)
    n = length(x)
    T = Matrix{ComplexF64}(undef, n, n)
    for i in 1:n, j in 1:n
        sum = 0.0 + 0.0im
        for (idx, m) in enumerate(fmodes)
            am = α + 2π * m / period
            sum += im * beta[idx] * exp(im * am * (x[i] - x[j]))
        end
        T[i, j] = sum / n
    end
    return T
end

"""
    apply_tbc!(Δ, sp::SamplingPoints, k; α = 0.0, homo = 1.0, period = 2π)

Impose the Transparent Boundary Condition for the top and bottom boundary.

# Arguments
- `Δ`: the matrix after imposing the periodic (quasi-periodic) boundary condition, see [`apply_bc`](@ref)
- `sp`: all sampling points
- `k`: the wavenumber
- `α`: the Bloch wavenumber
- `homo`: the dielectric constant of the homogeneous medium outside the periodic slab
- `period`: the period of the periodic slab
"""
function  apply_tbc!(Δ, sp::SamplingPoints, k; α = 0.0, homo = 1.0, period = 2π)
    n = sp.n
    # Get bottom and top x coordinates
    xb = @views sp.cc[1, 1:n]
    xt = @views sp.cc[1, 2n+1:3n]
    
    if iseven(n)
        fmodes = collect(-n÷2 : n÷2 - 1)
    else
        fmodes = collect(-(n-1)÷2 : (n-1)÷2)
    end
    
    β = [beta_m(m, k; α = α, homo = homo, period = period) for m in fmodes]
    
    Tb = assemble_tbc(xb, β, fmodes; α = α, period = period)
    Tt = assemble_tbc(xt, β, fmodes; α = α, period = period)
    
    @views Δ[1:n, 1:n] .-= Tb 
    @views Δ[2n+1:3n, 2n+1:3n] .-= Tt
    
    return Δ
end