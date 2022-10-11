using DrWatson
@quickactivate
using DynamicalSystems
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9



function compute_thomas(di::Dict)
    @unpack b, res = di
    ds = Systems.thomas_cyclical(b = b)
    xg = yg = zg = range(-7, 7,length=10000)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, safety_counter_max = Int(1e7), show_progress = true, diffeq)
    y1 = y2 = range(-5, 5, length = res)
    bsn = [ mapper([x,y,0.]) for x in y1, y in y2]
    att = mapper.bsn_nfo.attractors
    # bsn, att = basins_of_attraction(mapper, (y1,y2); show_progress = true)
    grid = (y1,y2)
    return @strdict(bsn, att, grid, b, res)
end


function print_fig(w, h, cmap, b, res)
    params = @strdict res b
    data, file = produce_or_load(
        datadir("basins"), params, compute_thomas;
        prefix = "thomas", storepatch = false, suffix = "jld2", force = false
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
    if isnothing(cmap)
        heatmap!(ax, xg, yg, bsn, rasterize = 1)
    else
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    end
    save(string(projectdir(), "/plots/thomas",res,".png"),fig)
end

b=0.1665
print_fig(600,600, nothing, b, 200) 
