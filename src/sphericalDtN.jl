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

"""
    BoundaryCache

# Fields
- `orders::Vector{Int64}`: the orders of Bessel functions and Neumann functions
- `J_ext::Matrix{Float64}`: Bessel functions at the sampling points of order in `orders`
- `Y_ext::Matrix{Float64}`: Neumann functions at the sampling points of order in `orders`
"""
struct BoundaryCache
    orders::Vector{Int64}
    J_ext::Matrix{Float64}
    Y_ext::Matrix{Float64}
end

"""
    build_boundary_cache(sp::SamplingPoints, Ne, k, ext = 1.0)

Build a `BoundaryCache` for computing the Dirichlet-to-Neuamnn map based on 
the spherical expansion.

# Arguments
- `sp`: the sampling points
- `Ne`: the number of terms in the truncated spherical expansion
- `k`: the wavenumber
- `ext`: the dielectric constant outside the cylinder

# Note 
- `J_ext` and `Y_ext` are matrices of Bessel functions and Neumann functions, respectively. Their rows are indexed by radii (the sampling points) and columns are indexed by orders.
- We uitilize the symmetry of the square. So this function is not applicable for general contours.
"""
function build_boundary_cache(sp::SamplingPoints, Ne, k, ext = 1.0)
    Ns = 4 * sp.n
    (Ns == Ne) || throw(ArgumentError("Ne must equal 4 * sp.n for a square matrix!"))

    Ne2 = Ne + 2
    k_ext = k * sqrt(ext)
    
    orders = Vector{Int64}(undef, Ne2)
    
    # Get unique radii due to symmetry
    Nu = cld(sp.n, 2) 
    ru = Vector{Float64}(undef, Nu)                                                                                                                                                                
    for i in 1:Nu                                                                                                                                                                                      
        ru[i] = sp.pc[1, i]
    end      
    
    # Compute Bessel functions and Neumann functions for unique radii
    Ju = Matrix{Float64}(undef, Nu, Ne2)
    Yu = similar(Ju)
    args = k_ext .* ru
    for j in 1:Ne2 
        orders[j] = j - 2 - (Ne ÷ 2)
        for i in 1:Nu
            Ju[i, j] = besselj(orders[j], args[i])
            Yu[i, j] = bessely(orders[j], args[i])
        end
    end

    # Compute the mapping between unique radii and all radii
    ri = Vector{Int64}(undef, Ns)
    for edge in 0:3
        offset = edge * sp.n
        for i in 1:sp.n
            ri[offset + i] = (i <= Nu ? i : (sp.n - i + 1))
        end
    end

    # Extend to all radii
    J_ext = Matrix{Float64}(undef, Ns, Ne2)
    Y_ext = similar(J_ext)
    for j = 1:Ne2
        for i = 1:Ns
            ii = ri[i]
            J_ext[i, j] = Ju[ii, j]
            Y_ext[i, j] = Yu[ii, j]
        end
    end

    return BoundaryCache(orders, J_ext, Y_ext)
end

"""
    assemble_dtn(sp::SamplingPoints, Ne, coeff, bc::BoundaryCache, k, ext = 1.0)

Compute the Dirichlet-to-Neumann matrix on a square based on the spherical expansion.

# Arguments
- `sp`: the sampling points
- `Ne`: the number of terms in the truncated spherical expansion
- `coeff`: from [`get_coeff()`](@ref), coefficients in front of Bessel functions and Neumann functions in the spherical expansion outside the cylinder
- `bc`: the `BoundaryCache`
- `k`: the wavenumber
- `ext`: the dielectric constant outside the cylinder
"""
function assemble_dtn(sp::SamplingPoints, Ne, coeff, bc::BoundaryCache, k, ext = 1.0)
    Ns = 4 * sp.n
    (Ns == Ne) || throw(ArgumentError("Ne must equal 4 * sp.n for a square matrix!"))
    (size(coeff, 2) == Ne) || throw(ArgumentError("coeff size mismatch"))

    k_ext = k * sqrt(ext)
    A = Matrix{ComplexF64}(undef, Ne, Ne)
    B = similar(A)
    
    for j = 1:Ne 
        order = bc.orders[j + 1]
        aj = coeff[1, j]
        bj = coeff[2, j]

        for i = 1:Ns
            # extract the polar coordinate
            r, θ = sp.pc[1, i], sp.pc[2, i]
            expo = exp(im * order * θ)
            
            Φ = aj * bc.J_ext[i, j + 1] + bj * bc.Y_ext[i, j + 1]
            
            # Derivatives
            Jp = (bc.J_ext[i, j] - bc.J_ext[i, j + 2]) / 2
            Yp = (bc.Y_ext[i, j] - bc.Y_ext[i, j + 2]) / 2
            Ψ = k_ext * (aj * Jp + bj * Yp)
            
            # extract the normal vector
            ν1, ν2 = sp.ν[1, i], sp.ν[2, i]
            c1 = cos(θ) * ν1 + sin(θ) * ν2
            c2 = cos(θ) * ν2 - sin(θ) * ν1
            
            A[i, j] = Φ * expo
            B[i, j] = Ψ * expo * c1 + Φ * (im * order * expo / r) * c2
        end
    end
    
    return B / A, A
end