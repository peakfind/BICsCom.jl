# empty declaration for function in extensions

"""
    plot_square!(ax, sq::Square)

Plot the square `sq`.

This function is provided by the `CMakieExt` extension and requires `CairoMakie`
to be loaded.
"""
function plot_square! end

"""
    plot_samplingpoints!(ax, sp::SamplingPoints;
                         skwargs = (marker = 'x', color = :red),
                         akwargs = (lengthscale=0.05, color=:red, shaftwidth=1, tipwidth=4))

Plot the sampling points `sp` with their corresponding normal vectors.

This function is provided by the `CMakieExt` extension and requires `CairoMakie`
to be loaded.

# Arguments
- `skwargs`: keyword arguments for `scatter!` in `CairoMakie`
- `akwargs`: keyword arguments for `arrows2d!` in `CairoMakie`
"""
function plot_samplingpoints! end

"""
    plot_min_svals!(ax, ks, msvals;
                    skwargs = (color = :tomato, marker = 'x', markersize = 5),
                    lkwargs = (color = :blue))

Plot the minimal singular values for frequencies in `ks`.

This function is provided by the `CMakieExt` extension and requires `CairoMakie`
to be loaded.

# Arguments
- `ks`: the range of frequencies we sweep
- `msvals`: the minimal singular values corresponding to the frequencies in `ks`
- `skwargs`: keyword arguments for `scatter!` in `CairoMakie`
- `lkwargs`: keyword arguments for `lines!` in `CairoMakie`
"""
function plot_min_svals! end

"""
    plot_field_mode!(ax, xs, ys, field; hkwargs = (colormap = :coolwarm))

Plot the field of the mode of a BIC.

# Arguments
- `xs`: discrete points along the ``x`` direction (the direction with periodicity)
- `ys`: discrete points along the ``y`` direction
- `field`: a matrix of size `length(xs)` times `length(ys)`
- `hkwargs`: keyword arguments for `heatmap!` in `CairoMakie`
"""
function plot_field_mode! end