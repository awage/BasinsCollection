using DrWatson
@quickactivate # exports Attractors, GLMakie and other goodies in `src`
using Attractors
using Attractors
using OrdinaryDiffEq
using CairoMakie
using LaTeXStrings


@inline @inbounds function duffing_pforced(u, p, t)
    a, b, c, ω = p; x,y = u
    du1 = y 
    du2 = a*x -b*x^2 -c*y + x*5*cos(ω*t)
    return SVector{2}(du1, du2)
end



function compute_basins_duffing_pforced(di::Dict)
    @unpack a, b, c, ω, res = di
    # diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    diffeq = (;reltol = 1e-9, alg = Tsit5(), maxiters = 1e6, dt = 0.01)
    ds = ContinuousDynamicalSystem(duffing_pforced, rand(2), [a,b,c, ω]; diffeq)
    xg = yg = range(-2.2,2.2,length = res)
    smap = StroboscopicMap(ds, 2*pi/ω)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


function print_fig(w,h, a, b, c, ω, res)
    params = @strdict a b c ω res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_duffing_pforced;
        prefix = "duffing_pforced", storepatch = false, suffix = "jld2", force = true
    )
    @unpack bsn, grid = data
    xg, yg = grid

    fig = Figure(size = (w, h))
    ax = Axis(fig[1,1], ylabel = L"$\dot{x}$", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, xg, yg, bsn, rasterize = 1)
    save(string("../plots/basins_duffing_pforced", res,".png"),fig)
end


a = b = 1; c = 0.2; ω=1.;  # smooth boundary
print_fig(600, 600, a, b, c, ω, 600) 
