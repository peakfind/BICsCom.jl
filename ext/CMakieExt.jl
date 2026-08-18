module CMakieExt

using BICsCom
using CairoMakie

export plot_square!, plot_samplingpoints!
export plot_min_svals!

"""
    plot_square!(ax, sq::Square)

Plot the square `sq`.
"""
function plot_square!(ax, sq::Square)
    r = sq.center[1] + sq.semi
    l = sq.center[1] - sq.semi
    t = sq.center[2] + sq.semi
    b = sq.center[2] - sq.semi
    lines!(ax, [r, l, l, r, r], [t, t, b, b, t])
    
    return ax
end


"""
    plot_samplingpoints!(ax, sp::SamplingPoints; 
                         skwargs = (marker = 'x', color = :red), 
                         akwargs = (lengthscale=0.05, color=:red, shaftwidth=1, tipwidth=4))

Plot the sampling points `sp` with corresponding normal vectors.

# Arguments
- `skwargs`: keyword arguments for `scatter` in CairoMakie
- `akwargs`: keywork arguments for `arrows2d` in CairoMakie
"""
function plot_samplingpoints!(ax, sp::SamplingPoints; skwargs = (marker = 'x', color = :red), akwargs = (lengthscale=0.05, color=:red, shaftwidth=1, tipwidth=4))
    scatter!(ax, sp.cc[1,:], sp.cc[2,:], skwargs...)
    arrows2d!(ax, sp.cc[1,:], sp.cc[2,:], sp.ν[1,:], sp.ν[2,:], akwargs...)
    
    return ax
end

"""
    plot_min_svals!(ax, ks, msvals; 
                    skwargs = (color = :tomato, marker = 'x', markersize = 5), 
                    lkwargs = (color = :blue))

Plot the minimal singular values for frequencies in `ks`.

# Arguments
- `ks`: the range of frequencies we sweep
- `msvals`: the minimal singular values corresponding to the frequencies in `ks`
- `skwargs`: keyword arguments for `scatter` in CairoMakie
- `lkwargs`: keyword arguments for `lines` in CairoMakie
"""
function plot_min_svals!(ax, ks, msvals; skwargs = (color = :tomato, marker = 'x', markersize = 5), lkwargs = (color = :blue))
    scatter!(ax, ks, msvals, skwargs...)
    lines!(ax, ks, msvals, lkwargs...)
    
    return ax
end

end