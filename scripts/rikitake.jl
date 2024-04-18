using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter
# This system is problematic: very very long transient (t > 2000 sometimes) that can be mistaken with attractors.




function compute_rikitake(di::Dict)
    @unpack μ, α, res = di
    ds = Systems.rikitake(μ = μ, α = α)
    xg = yg = zg = range(-5,5,length=10000)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    # pmap = poincaremap(ds, (3, 0.); rootkw = (xrtol = 1e-8, atol = 1e-8), diffeq)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, safety_counter_max = Int(1e7), show_progress = true)
    y1 = y2 = range(-2.5, 2.5, length = res)
    bsn = @showprogress [ mapper([x,y,0.]) for x in y1, y in y2]
    att = mapper.bsn_nfo.attractors
    # bsn, att = basins_of_attraction(mapper, (y1,y2); show_progress = true)
    grid = (y1,y2)
    return @strdict(bsn, att, grid, μ, α, res)
end


function print_fig(w, h, cmap, μ, α, res)
    params = @strdict res μ α
    data, file = produce_or_load(
        datadir("basins"), params, compute_rikitake;
        prefix = "rikitake", storepatch = false, suffix = "jld2", force = false
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
    if isnothing(cmap)
        heatmap!(ax, xg, yg, bsn, rasterize = 1)
    else
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    end
    save(string(projectdir(), "/plots/rikitake",res,".png"),fig)
end

print_fig(600,600, nothing, 0.47, 1., 700) 
