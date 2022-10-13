using DrWatson
@quickactivate "PerspectiveFigures" # exports DynamicalSystems, GLMakie and other goodies in `src`
using DynamicalSystems
using OrdinaryDiffEq:Vern9
using CairoMakie
using LaTeXStrings


@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end



function compute_basins_duffing(di::Dict)
    @unpack d, F, ω, res = di

    ds = ContinuousDynamicalSystem(duffing, rand(2), [d, F, ω])
    xg = yg = range(-2.2,2.2,length = res)
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    smap = stroboscopicmap(ds, 2*pi/ω; diffeq)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    bsn, att = basins_of_attraction(mapper)
    grid = (xg,yg)
    return @strdict(bsn, att, grid, d, F, ω, res)
end


function print_fig(w,h,cmap, d, F, ω, res)
    params = @strdict d F ω res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_duffing;
        prefix = "duffing", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    xg, yg = grid

    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"$\dot{x}$", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    save(string("../plots/basins_duffing", d,".png"),fig)

end

function get_Sb(d, F, ω, res)
    params = @strdict d F ω res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_duffing;
        prefix = "duffing", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    basins_fractal_test(bsn; ε = 10)
    return basin_entropy(bsn)
end

cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
d = 0.1; F=0.1; ω=0.1;  # smooth boundary
print_fig(600, 500, cmap, d, F, ω, 600) 

d = 0.4; F=0.1; ω=0.1;  # smooth boundary
print_fig(600, 500, cmap, d, F, ω, 600) 
