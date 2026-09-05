# # Standing waves for dielectric cylinders under TE polarization

# ## Parameters
# This example is from [Hu2015](@cite)

# First we need to import needed packages
using BICsCom
using CairoMakie
using LinearAlgebra

# Then we specify the parameters 
n = 5
inn = 11.6
ext = 1.0
hom = 1.0
r = 0.3 * 2π

# We construct a empty for minimal singular values 
msv = Float64[]

# a square with period ``2\pi``
sq = Square([0.0, 0.0], π)

# generate sampling points along the boundary of the square 
sp = get_samplingpoints(sq, n)

# We sweep the frequency 
ks = 0.4:1e-5:1.0

for k in ks
    cydc = build_cylinder_cache(4n, k, r, inn, ext)
    coef = get_coeff(4n, cydc, inn, ext)
    bydc = build_boundary_cache(sp, 4n, k, ext)
    Λ, _ = assemble_dtn(sp, 4n, coef, bydc, k, ext)
    Δ = apply_bc(Λ, n)
    apply_tbc!(Δ, sp, k; homo = hom)
    
    s = svdvals(Δ)
    push!(msv, minimum(s))
end

# We plot all minimal singular values for all frequencies in `ks`
with_theme(theme_latexfonts()) do
    fig = Figure()
    axi = Axis(fig[1, 1], yscale = log10, xminorgridvisible = true, yminorgridvisible = true, xminorgridstyle = :dash, 
               yminorgridstyle = :dash, xminorticks = IntervalsBetween(10), yminorticks = IntervalsBetween(20)) 
    plot_min_svals!(axi, ks, msv)
    
    fig
end