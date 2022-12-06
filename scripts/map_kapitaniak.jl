using DrWatson
@quickactivate
using DynamicalSystems
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9

# Bifurcations from locally to globally riddled basins
# T. Kapitaniak, Yu. Maistrenko, A. Stefanski, and J. Brindley
# https://doi.org/10.1103/PhysRevE.57.R6253
function kapitaniak_map!(dz, z, p, n)
    xn = z[1]; yn = z[2]
    l,pp,d1,d2 = p
    dz[1] = pp*xn + l/2*(1-pp/l)*(abs(xn+1/l)-abs(xn-1/l))+d1*(yn-xn)
    dz[2] = pp*yn + l/2*(1-pp/l)*(abs(yn+1/l)-abs(yn-1/l))+d2*(xn-yn)
    return
end

# dummy function to keep the initializator happy
function kapitaniak_map_J(J, z0, p, n)
    return
end



function compute_kapitaniak(di::Dict)
    @unpack l, pp, d1, d2, res = di
    ds = DiscreteDynamicalSystem(kapitaniak_map!, [1.0, 0.0], [l, pp, d1, d2], kapitaniak_map_J)
    yg = xg = range(-3., 3., length = 10000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, safety_counter_max = Int(1e7), show_progress = true)
    yg = xg = range(-2, 2, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end


function print_fig(w, h, cmap, l, pp, d1, d2, res)
    params = @strdict res l pp d1 d2
    data, file = produce_or_load(
        datadir("basins"), params, compute_kapitaniak;
        prefix = "kapitaniak", storepatch = false, suffix = "jld2", force = true
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
    save(string(projectdir(), "/plots/kapitaniak",res,".png"),fig)
end

# l = 1.3; pp = -2.; d1 = 0.725; d2 = 0.725;
l = √2; pp = -√2; d1 = d2 = -0.935;

print_fig(600,600, nothing, l, pp, d1, d2, 2000) 
