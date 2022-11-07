using DrWatson
@quickactivate 
using OrdinaryDiffEq:Vern9
using DynamicalSystems
using CairoMakie
using ProgressMeter
# Sommerer, J. C. (1995). The end of classical determinism. Johns Hopkins APL Technical Digest, 16(4), 333.
function forced_particle!(du, u, p, t)
    γ = 0.632; f₀ = 1.0688  ; ω = 2.2136; 
    s = 20.; p = 0.098; k = 10.;
    x₀=1. ; y₀=0.;
    x, y, dx, dy = u
    du[1] = dx
    du[2] = dy
    du[3] = -γ*dx -(-4x*(1-x^2) + 2*s*x*y^2) +  f₀*sin(ω*t)*x₀
    du[4] = -γ*dy -(2*y*s*(x^2-p)+4*k*y^3) +  f₀*sin(ω*t)*y₀
end


function _get_basins_sommerer(d)
    @unpack res = d
    xg = range(-3,3,length = 3000)
    yg = range(-3,3,length = 3000)
    df = ContinuousDynamicalSystem(forced_particle!,rand(4),(0.0,20.0))
    diffeq = (reltol = 1e-9,  alg = Vern9())
    ω = 2.2136
    smap = stroboscopicmap(df, 2π/ω; diffeq)
    # psys = projected_integrator(smap, [1,2], [0., 0,])
    mapper = AttractorsViaRecurrences(smap, (xg, yg, xg, yg); 
            # mx_chk_att = 5,
            sparse = true)

    xg = range(-1.,1.,length=res)
    yg = range(-1.,1.,length=res)
    basins = @showprogress [ mapper([x,y,0., 0.]) for x in xg, y in yg]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(basins, xg, yg, att)
end



function print_fig(w,h,cmap, res)

    # res = 1500
    data, file = produce_or_load(
        datadir("basins"), # path
        @dict(res), # container
        _get_basins_sommerer, # function
        prefix = "basin_sommerer", # prefix for savename
        force = true
    )
    @unpack basins, xg, yg = data

    # xg = range(-2,2,length = res)
    # yg = range(0.,2,length = res)
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"y_0", xlabel = L"x_0", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    if isnothing(cmap)
        heatmap!(ax, xg, yg, basins, rasterize = 1)
    else
        heatmap!(ax, xg, yg, basins, rasterize = 1, colormap = cmap)
    end
    save("../plots/basins_riddle_sommerer.png",fig)
end



# cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
print_fig(600, 600, nothing, 1000) 
