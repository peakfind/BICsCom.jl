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
    axi = Axis(fig[1, 1], yscale = log10, 
               xminorgridvisible = true, yminorgridvisible = true, 
               xminorgridstyle = :dash, yminorgridstyle = :dash, 
               xminorticks = IntervalsBetween(10), 
               yminorticks = IntervalsBetween(20), 
               title = "minimal singular values vs. frequencies") 
    plot_min_svals!(axi, ks, msv)
    fig
end

# Next we plot the fields of modes associated to BICs we computed above.
# For the BIC's frequency ``k1 = 0.4112``, 
k1 = 0.4112
mf1 = compute_mode(sp, k1, r, inn, ext, hom)
# We compute the field on a grid `xs` times `ys`
xs = -π:0.05:π
ys = -π:0.05:π
f1 = evaluate_field(mf1, xs, ys)

with_theme(theme_latexfonts()) do
    fig2 = Figure()
    
    angles = range(0, 2π, length=200)
    xc = r .* cos.(angles)
    yc = r .* sin.(angles)
    
    axi2 = Axis(fig2[1, 1], aspect = 1, xlabel = L"$x$", ylabel = L"$y$")
    hm, _ = plot_field_mode!(axi2, xs, ys, real.(f1))
    lines!(axi2, xc, yc, color = :black)
    Colorbar(fig2[1, 2], hm)
    fig2
end