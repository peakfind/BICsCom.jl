# # Propagating waves for dielectric cylinders under TE polarization

# ## Problem 
#
# In [Yuan2017](@cite), the authors report two propagating BICs:
# - x-even BIC with the Bloch wave number ``\beta L/ 2\pi = 0.0776`` and the frequency ``\omega L / 2\pi c = 0.4854``
# - x-odd BIC with the Bloch wave number ``\beta L/ 2\pi = 0.2483`` and the frequency ``\omega L / 2\pi c = 0.6702``

using BICsCom
using CairoMakie
using LinearAlgebra

# We first specify the parameters we need
n = 10
inn = 11.56
ext = 1.0
hom = 1.0
r = 0.35 * 2π
α = 0.2483

# collect minimal singular values
msv = Float64[]

# a square with period 2π
sq = Square([0.0, 0.0], π)

# generate sampling points along the boundary of the square
sp = get_samplingpoints(sq, n)

# sweep the frequency
ks = 0.4:1e-5:0.8

for k in ks
    cydc = build_cylinder_cache(4n, k, r, inn, ext)
    coef = get_coeff(4n, cydc, inn, ext)
    bydc = build_boundary_cache(sp, 4n, k, ext)
    Λ, _ = assemble_dtn(sp, 4n, coef, bydc, k, ext)
    Δ = apply_bc(Λ, n; α = α)
    apply_tbc!(Δ, sp, k; α = α, homo = hom)
    
    s = svdvals(Δ)
    push!(msv, minimum(s))
end

# Plot 
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

# We extract the frequency for this x-odd BIC
ind = argmin(msv)
ko = ks[ind]

mfo = compute_mode(sp, ko, r, inn, ext, hom; α = α, period = 2π)

xs = -π:0.05:π
ys = -π:0.05:π
fo = evaluate_field(mfo, xs, ys)

with_theme(theme_latexfonts()) do 
    fig2 = Figure()
    
    angles = range(0, 2π, length=200)
    xc = r .* cos.(angles)
    yc = r .* sin.(angles)
    
    axi2 = Axis(fig2[1, 1], aspect = 1, xlabel = L"$x$", ylabel = L"$y$")
    hm, _ = plot_field_mode!(axi2, xs, ys, real.(fo))
    lines!(axi2, xc, yc, color = :black)
    Colorbar(fig2[1, 2], hm)
    fig2
end