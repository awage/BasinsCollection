using DrWatson
@quickactivate
using DynamicalSystems
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9

# International Journal of Bifurcation and ChaosVol. 03, No. 04, pp. 803-842 (1993) Tutorials and ReviewsNo Access
# THE BOGDANOV MAP: BIFURCATIONS, MODE LOCKING, AND CHAOS IN A DISSIPATIVE SYSTEM
# DAVID K. ARROWSMITH, JULYAN H. E. CARTWRIGHT, ALEXIS N. LANSBURY, and COLIN M. PLACE
# https://doi.org/10.1142/S021812749300074X
function bogdanov_map!(dz, z, p, n)
    x, y = z
    μ, k, ε = p
    dz[2] = y + ε*y + k*x*(x-1) + μ*x*y
    dz[1] = x + dz[2]
    return
end

# dummy function to keep the initializator happy
function bogdanov_map_J(J, z0, p, n)
    return
end



function compute_bogdanov(di::Dict)
    @unpack μ, k, ε, res = di
    ds = DiscreteDynamicalSystem(bogdanov_map!, [1.0, 0.0], [μ, k, ε], bogdanov_map_J)
    yg = xg = range(-2., 2., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        mx_chk_fnd_att = 1000,
        mx_chk_loc_att = 1000, mx_chk_safety = Int(1e7), show_progress = true)
    yg = xg = range(-0.8, 1.1, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end


function print_fig(w, h, cmap, μ, k, ε, res)
    params = @strdict res μ k ε
    data, file = produce_or_load(
        datadir("basins"), params, compute_bogdanov;
        prefix = "bogdanov", storepatch = false, suffix = "jld2", force = true
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
    save(string(projectdir(), "/plots/bogdanov",res,".png"),fig)
end

μ = -0.1
k = 1.2
ε = 0.0125 
print_fig(600,600, nothing, μ, k, ε, 1000) 
