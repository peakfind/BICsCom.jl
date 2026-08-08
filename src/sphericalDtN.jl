"""
    CylinderCache
    
Cache for Bessel functions and Neumann function required in the spherical 
expansion outside the cylinder.
    
# Fields
- `r`: the radius of the cylinder
- `orders`: the orders of Bessel functions and Neumann functions
- `J_inn`: Bessel functions inside the cylinder of order in `orders`
- `J_ext`: Bessel functions outside the cylinder of order in `orders`
- `Y_ext`: Neumann functions outside the cylinder of order in `orders`
"""
struct CylinderCache
    r::Float64
    orders::Vector{Int64}
    J_inn::Vector{Float64}
    J_ext::Vector{Float64}
    Y_ext::Vector{Float64}
end

"""
    build_cylinder_cache(Ne, k, r, inn = 1.0, ext = 1.0)

Build a `CylinderCache` for computing the coefficients in the spherical expansion 
outside the cylinder.

# Arguments
- `Ne`: the number of terms in the truncated spherical expansion
- `k`: the wavenumber
- `r`: the radius of the cylinder
- `inn`: the dielectric constant inside the cylinder
- `ext`: the dielectric constant outside the cylinder
"""
function build_cylinder_cache(Ne, k, r, inn = 1.0, ext = 1.0) 
    Ne2 = Ne + 2
    t_inn = k * sqrt(inn) * r
    t_ext = k * sqrt(ext) * r
    
    orders = Vector{Int64}(undef, Ne2)
    J_inn = Vector{Float64}(undef, Ne2)
    J_ext = similar(J_inn)
    Y_ext = similar(J_inn)
    
    for i = 1:Ne2
        o = i - 2 - (Ne ÷ 2)
        orders[i] = o
        J_inn[i] = besselj(o, t_inn) 
        J_ext[i] = besselj(o, t_ext)
        Y_ext[i] = bessely(o, t_ext)
    end
    
    return CylinderCache(r, orders, J_inn, J_ext, Y_ext)
end

"""
    get_coeff(Ne, cydc::CylinderCache, inn = 1.0, ext = 1.0)

Compute the coefficients in the spherical expansion outside the cylinder.

# Arguments
- `Ne`: the number of terms in the truncated spherical expansion
- `cydc`: the `CylinderCache`
- `inn`: the dielectric constant inside the cylinder
- `ext`: the dielectric constant outside the cylinder
"""
function get_coeff(Ne, cydc::CylinderCache, inn = 1.0, ext = 1.0)
    sqr_inn = sqrt(inn)
    sqr_ext = sqrt(ext)

    Ji = cydc.J_inn
    Je = cydc.J_ext
    Ye = cydc.Y_ext
    
    coeff = Matrix{Float64}(undef, 2, Ne)
    
    for i in 1:Ne
        lhs = [Je[i + 1] Ye[i + 1]; (Je[i] - Je[i + 2]) (Ye[i] - Ye[i + 2])]
        rhs = [Ji[i + 1], sqr_inn * (Ji[i] - Ji[i + 2]) / sqr_ext]
        coeff[:, i] = lhs \ rhs
    end

    return coeff
end