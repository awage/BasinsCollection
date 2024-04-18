using DrWatson
@quickactivate # exports Attractors, GLMakie and other goodies in `src`
using Attractors
using Attractors
using OrdinaryDiffEq
using CairoMakie
using LaTeXStrings

# φ = (- b  ̇φ - c φ - m*L*(μ − A_2*ω^2 cos ωt) sin φ - m*L*A_1*ω^2 sin ωt cos φ)/(m*L^2) 
# Parametrically excited pendulum systems with several equilibrium positions. bifurcation analysis and rare attractors
# Alex v. Klokov∗ and Mikhail v. Zakrzhevsky
# http://dx.doi.org/10.1142/S0218127411030167
@inline @inbounds function pend_pforced(u, p, t)
    A_1, A_2, ω = p; x,y = u
    μ = 9.81; L = m = c = 1; b = 0.2
    du1 = y 
    # du2 = (-b*y -c*x -m*L*(μ - A_2*ω^2*cos(ω*t))*sin(x) - m*L*A_1*ω^2*sin(ω*t)*cos(x))/(m*L^2)
    du2 = -b*y -x -(μ - A_2*ω^2*cos(ω*t))*sin(x) - A_1*ω^2*sin(ω*t)*cos(x)
    return SVector{2}(du1, du2)
end



function compute_basins_pend_pforced(di::Dict)
    @unpack A_1, A_2, ω, res = di

    ds = ContinuousDynamicalSystem(pend_pforced, rand(2), [A_1, A_2, ω])
    xg = yg = range(-14,14,length = 4001)
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e9)
    smap = stroboscopicmap(ds, 2*pi/ω; diffeq)
    mapper = AttractorsViaRecurrences(smap, (xg, yg);
        mx_chk_fnd_att = 100, mx_chk_loc_att = 100, mx_chk_safety = Int(1e9), show_progress = true)
    xg = yg = range(-10,10,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    grid = (xg,yg)
    return @strdict(bsn, att, grid, res)
end


function print_fig(w,h, A_1, A_2, ω, res)
    params = @strdict A_1 A_2 ω res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_pend_pforced;
        prefix = "pend_pforced", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    xg, yg = grid

    fig = Figure(size = (w, h))
    ax = Axis(fig[1,1], ylabel = L"$\dot{\varphi}$", xlabel = L"\varphi", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, xg, yg, bsn, rasterize = 1)
    save(string("../plots/basins_pend_pforced", res,".png"),fig)
end

# A_1 = 0.5; A_2 = 0; ω = 1.5
A_1 = 0.; A_2 = 4.1; ω = 1.5
print_fig(600, 600, A_1, A_2, ω, 500) 
