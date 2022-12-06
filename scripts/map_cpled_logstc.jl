using DrWatson
@quickactivate
using DynamicalSystems
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9

# Transverse instability and riddled basins in a system of two coupled logistic maps
# Yu. L. Maistrenko, V. L. Maistrenko, and A. Popovich E. Mosekilde
#https://doi.org/10.1103/PhysRevE.57.2713
function cplog_map(dz, z, p, n)
    xn = z[1]; yn = z[2]
    a, ε = p
    dz[1] = a*xn*(1-xn) + ε*(yn-xn) 
    dz[2] =  a*yn*(1-yn) + ε*(xn-yn)   
    return
end

# dummy function to keep the initializator happy
function cplog_map_J(J, z0, p, n)
    return
end



function compute_cplog(di::Dict)
    @unpack a, ε, res = di
    ds = DiscreteDynamicalSystem(cplog_map, [1.0, 0.0], [a, ε], cplog_map_J)
    yg = xg = range(-2., 2., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, safety_counter_max = Int(1e7), show_progress = true)
    yg = xg = range(-0.05, 1.3, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, μ, j, res)
end


function print_fig(w, h, cmap, a, ε, res)
    params = @strdict res a ε
    data, file = produce_or_load(
        datadir("basins"), params, compute_cplog;
        prefix = "cplog", storepatch = false, suffix = "jld2", force = true
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
    save(string(projectdir(), "/plots/cplog",res,".png"),fig)
end

# a = 3.57480493875920; ε = -0.2
a = 3.6; ε = -1.
print_fig(600,600, nothing, a, ε, 1000) 
