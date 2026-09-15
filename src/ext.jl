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
    plot_field_mode!(ax, xs, ys, field; transform = abs, hkwargs = (; colormap = :coolwarm))
    plot_field_mode!(ax, mf::ModeField, sq::Square, xs, ys; transform = abs, gauge = true, hkwargs = (; colormap = :coolwarm))

Plot the field of the mode of a BIC.

The first method plots the matrix `field` of size `length(xs)` times `length(ys)`.
The second method evaluates the field of `mf` on the grid with Bloch continuation,
see [`evaluate_field`](@ref).

# Keyword Arguments
- `transform = abs`: the function applied elementwise to the complex field, e.g.
  `abs`, `real`, or `imag`
- `gauge = true`: apply [`gauge_field!`](@ref) before `transform`. With the
  default `transform = abs` this is a visual no-op, and it makes
  `transform = real` show the gauge-corrected standing-wave pattern for
  standing BICs (`α = 0`)
- `hkwargs`: keyword arguments for `heatmap!` in `CairoMakie`

This function is provided by the `CMakieExt` extension and requires `CairoMakie`
to be loaded.
"""
function plot_field_mode! end