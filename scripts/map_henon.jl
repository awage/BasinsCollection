using DrWatson
@quickactivate
using DynamicalSystems
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9



function henon_map(dz, z, p, n)
    xn = z[1]; yn = z[2]
    μ, j = p
    dz[1] = 1 -μ*xn^2 +yn
    dz[2] = -j*xn  
    return
end

# dummy function to keep the initializator happy
function henon_map_J(J, z0, p, n)
    return
end



function compute_henon(di::Dict)
    @unpack μ, j, res = di
    ds = DiscreteDynamicalSystem(henon_map, [1.0, 0.0], [μ, j], henon_map_J)
    yg = xg = range(-2., 2., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, safety_counter_max = Int(1e7), show_progress = true)
    yg = xg = range(-2, 2, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, μ, j, res)
end


function print_fig(w, h, cmap, μ, j, res)
    params = @strdict res μ j
    data, file = produce_or_load(
        datadir("basins"), params, compute_henon;
        prefix = "henon", storepatch = false, suffix = "jld2", force = false
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
    save(string(projectdir(), "/plots/henon",res,".png"),fig)
end

μ = 1.08
j = 0.9
print_fig(600,600, nothing, μ, j, 1000) 
