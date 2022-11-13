using DrWatson
@quickactivate
using DynamicalSystems
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter


# International Journal of Bifurcation and Chaos, Vol. 26, No. 14 (2016) 1650233 (11 pages)
#DOI: 10.1142/S0218127416502333
#Crisis in Amplitude Control Hides in Multistability
#Int. J. Bifurcation Chaos 2016.26.
# Chunbiao Li, Julien Clinton Sprott, Hongyan Xing. 
function li_sprott!(du, u, p, t)
    @inbounds begin
    a = p[1]; b = p[2];
    x, y, z = u
    du[1] = y + y*z
    du[2] = y*z - a*x*z
    du[3] = b*z^2 - y^2
    end
end


function compute_li_sprott(di::Dict)
    @unpack a, b, res = di
    ds = ContinuousDynamicalSystem(li_sprott!, rand(3), [a,b])
    xg = yg = zg = range(-5,5,length=10000)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); sparse = true,    
        mx_chk_fnd_att = 1000,
        mx_chk_loc_att = 1000, safety_counter_max = Int(1e8), show_progress = true)
    x1 = range(-1, 1, length = res) 
    y1 = range(-5, 5, length = res)
    bsn = @showprogress [ mapper([x,y,-1.]) for x in x1, y in y1]
    att = mapper.bsn_nfo.attractors
    # bsn, att = basins_of_attraction(mapper, (y1,y2); show_progress = true)
    grid = (x1,y1)
    return @strdict(bsn, att, grid, res)
end


function print_fig(w, h, cmap, a, b, res)
    params = @strdict res a b
    data, file = produce_or_load(
        datadir("basins"), params, compute_li_sprott;
        prefix = "li_sprott", storepatch = false, suffix = "jld2", force = true
    )
    @unpack bsn, grid = data
    xg, yg = grid
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"$y$", xlabel = L"x", yticklabelsize = 30, 
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
    save(string(projectdir(), "/plots/li_sprott",res,".png"),fig)
end

print_fig(600,600, nothing, 13., 0.55, 500) 
