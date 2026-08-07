"""
    Square
    
A square in ``\\mathbb{R}^2``. The boundary of the square is the boundary of 
a periodic cell. We establish a Dirichlet-to-Neumann (DtN) map on the boundary based 
on the spherical expansion of a solution to the Helmholtz equation.
    
# Fields
- `center::Vector{Float64}`: the center of the square
- `semi::Float64`: the semi length of the edge
"""
struct Square
    center::Vector{Float64}
    semi::Float64
end

"""
    SamplingPoints

The sampling points along the boundary of a square. The points are placed along the 
boundary counter clockwisely (bottom --> left --> top --> right).

# Fields
- `n`: the number of the sampling points on one edge
- `cc`: the Cartesian coordinates of all sampling points
- `pc`: the Polar coordinates of all sampling points
- `ν`: the outer normal derivative of all sampling points
"""
struct SamplingPoints
    n::Int64
    cc::Matrix{Float64}
    pc::Matrix{Float64}
    ν::Matrix{Float64}
end

"""
    cartesian2polar(cc::Matrix{Float64})
    cartesian2polar!(pc::Matrix{Float64}, cc::Matrix{Float64})

Change coordinates from Cartesian to Polar.
"""
function cartesian2polar(cc::Matrix{Float64})
    N = size(cc, 2)
    pc = similar(cc)
    
    for i in 1:N 
        x, y = cc[1, i], cc[2, i]
        pc[1, i] = sqrt(x^2 + y^2)
        pc[2, i] = atan(y, x)
    end
    
    return pc
end

function cartesian2polar!(pc::Matrix{Float64}, cc::Matrix{Float64})
    N = size(cc, 2)
    
    for i in 1:N 
        x, y = cc[1, i], cc[2, i]
        pc[1, i] = sqrt(x^2 + y^2)
        pc[2, i] = atan(y, x)
    end
    
    return pc
end

"""
    get_samplingpoints(sq::Square, n::Int64)

Generate the sampling points along the bounary of a square `sq`.

# Arguments
- `sq::Square`: a square
- `n::Int64`: the number of sampling points on each edge
"""
function get_samplingpoints(sq::Square, n::Int64)
    # Precompute common constants
    N = 4n
    semi = sq.semi
    δ = (2 * semi) / n
    start = semi - δ/2

    cc = Matrix{Float64}(undef, 2, N)
    pc = similar(cc)
    ν = similar(cc)
    
    # Construct the sampling points in Cartesian Coordinate
    # bottom edge (x: left <-- right, y = -semi) 
    for i = 1:n 
        cc[1,i] = start - (i - 1) * δ 
        cc[2,i] = -semi
    end

    # left edge (x = -semi, y: bottom --> top)
    for i = n+1:2n 
        cc[1,i] = -semi
        cc[2,i] = -start + (i - n - 1) * δ
    end

    # top edge (x: left --> right, y = semi)
    for i = 2n+1:3n 
        cc[1,i] = -start + (i - 2n - 1) * δ
        cc[2,i] = semi
    end

    # right edge (x = semi, y: bottom <-- top)
    for i = 3n+1:N
        cc[1,i] = semi
        cc[2,i] = start - (i - 3n - 1) * δ
    end
    
    # Construct pc by converting cc into Polar Coordinate
    cartesian2polar!(pc, cc)
    
    # Construct ν
    for i = 1:n 
        ν[1, i] = 0.0
        ν[2, i] = -1.0
    end
    
    for i = n+1:2n 
        ν[1, i] = -1.0
        ν[2, i] = 0.0
    end
    
    for i = 2n+1:3n 
        ν[1, i] = 0.0
        ν[2, i] = 1.0
    end
    
    for i = 3n+1:N
        ν[1, i] = 1.0
        ν[2, i] = 0.0
    end
    
    return SamplingPoints(n, cc, pc, ν)
end