# # Standing waves for slabs with air-holes under TE polarization
# This example comes from Fig. 2 in [Zhang2024](@cite)
using BICsCom
using CairoMakie
using LinearAlgebra

n = 9
inn = 1.0
ext = 8.2
hom = 1.0
r = 0.42 * 2π

msv = Float64[]

# period = 2π 
sq = Square([0.0, 0.0], π)
sp = get_samplingpoints(sq, n)

ks =  0.7:1e-5:0.9

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