using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors
using ChaosTools
using StaticArrays

# https://doi.org/10.48550/arXiv.2211.06921
# Disipative nontwist system
function dsnm!(dz, z, p, n)
    x,y = z
    γ = 0.1; 
    a, b = p
    dz[2] = (1-γ)*y - b*sin(2π*x)
    dz[1] = rem(x + a*(1-dz[2]^2),1, RoundNearest)
end



function compute_dsnm(di)
    @unpack a, b, res = di
    u0 = [0.; 0.]
    p = [a, b]
    df = DiscreteDynamicalSystem(dsnm!, u0, p) 
    x1 = range(-2, 2, length = 10001)
    y1 = range(-5, 5, length = 10001)
    grid_rec = (x1, y1)
    mapper = AttractorsViaRecurrences(df, grid_rec,
            # mx_chk_lost = 10, 
            mx_chk_fnd_att = 1000, 
            mx_chk_loc_att = 1000, 
            mx_chk_att = 4,
            safety_counter_max = Int(1e8),
            sparse = true, Ttr = 300)
    x = range(0., 1, length = res)
    y = range(-2, 2, length = res)
    grid = (x,y)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, res)
end



function print_fig(w, h, cmap, a, b, res) 
    params = @strdict a b res

    data, file = produce_or_load(
        datadir("basins"), params, compute_dsnm;
        prefix = "dsnm", storepatch = false, suffix = "jld2", force = true
    )


    @unpack bsn, grid = data
    xg ,yg = grid
    fig = Figure(size = (w, h))
    ax = Axis(fig[1,1], ylabel = L"y_0", xlabel = L"x_0", yticklabelsize = 30, 
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
    save(string("../plots/dsnm", res, ".png"),fig)
end


 
a = 0.55; b = 0.45 
print_fig(600,600, nothing, a, b, 1000) 

# p = [m, k]
# df = DiscreteDynamicalSystem(predator_prey!, rand(2), p) 
# u = trajectory(df, 1000, [0.3, 0.2])

