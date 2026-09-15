# # Standing waves for slabs with air-holes under TM polarization

# This example is from [Hu2015](@cite)

using BICsCom
using CairoMakie
using LinearAlgebra

n = 9
inn = 1.0
ext = 11.6
hom = 1.0
r = 0.3 * 2π

msv = Float64[];

sq = Square([0.0, 0.0], π) # period = 2π 
sp = get_samplingpoints(sq, n) # we have 4n points on the boundary

ks =  0.75:1e-5:0.8;

# We note that we need to specify `mode = :tm` in [`assemble_dtn`](@ref).
for k in ks
    cydc = build_cylinder_cache(4n, k, r, inn, ext)
    coef = get_coeff(4n, cydc, inn, ext, mode = :tm)
    bydc = build_boundary_cache(sp, 4n, k, ext)
    Λ, _ = assemble_dtn(sp, 4n, coef, bydc, k, ext, mode = :tm, homo = 1.0)
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

# We extract the BIC with the minimal frequency
ind = argmin(msv)
k1 = ks[ind]

mf1 = compute_mode(sp, k1, r, inn, ext, hom; mode = :tm)

xs = -π:0.05:π
ys = -4:0.05:4
fo = evaluate_field(mf1, xs, ys)

with_theme(theme_latexfonts()) do 
    fig2 = Figure()
    
    angles = range(0, 2π, length=200)
    xc = r .* cos.(angles)
    yc = r .* sin.(angles)
    
    axi2 = Axis(fig2[1, 1], xlabel = L"$x$", ylabel = L"$y$", aspect = DataAspect())
    hm, _ = plot_field_mode!(axi2, xs, ys, real.(fo))
    lines!(axi2, xc, yc, color = :black)
    lines!(axi2, [-π, π], [π, π], color = :black)
    lines!(axi2, [-π, π], [-π, -π], color = :black)
    Colorbar(fig2[1, 2], hm)
    fig2
end