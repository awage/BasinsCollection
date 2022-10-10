using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using DynamicalSystems

# Reference: Map with more than 100 coexisting low-period periodic attractors, Ulrike Feudel,  Celso Grebogi, Brian R. Hunt, and James A. Yorke
# PHYSICAL REVIEW E, VOLUME 54, NUMBER 1 1996
function kicked_feudel!(dz, z, p, n)
    xn = z[1]; yn = z[2]
    ν = p[1]; f0 = p[2]
    dz[1] = mod2pi(xn + yn)
    dz[2] = (1 - ν)*yn + f0*sin(xn + yn)
    return
end


function compute_kicked_rotor(di)
    @unpack f0, ν, res = di
    u0 = [0., 0.6]
    df = DiscreteDynamicalSystem(kicked_feudel!, u0, [ν,f0]) 
    xg = range(0, 2π, length = 10000)
    yg = range(-220, 220, length = 10000)
    grid_rec = (xg, yg)
    mapper = AttractorsViaRecurrences(df, grid_rec,
            mx_chk_lost = 10, 
            mx_chk_fnd_att = 3000, 
            mx_chk_loc_att = 3000, 
            mx_chk_att = 5,
            sparse = true
            )
    xg = range(0, 2π, length = res)
    yg = range(-π, π, length = res)
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, ν, f0 , res)
end



function print_fig(w, h, cmap, f0, ν, res) 
    params = @strdict f0 ν res

    data, file = produce_or_load(
        datadir("basins"), params, compute_kicked_rotor;
        prefix = "kicked_rotor", storepatch = false, suffix = "jld2", force = false
    )


    @unpack bsn, grid = data
    xg ,yg = grid
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"\dot{\theta}", xlabel = L"\theta", yticklabelsize = 30, 
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
    save(string("../plots/kicked_rotor_", res, ".svg"),fig)
end



res = 500
f0 = 4.
ν = 0.02
print_fig(600,600, nothing, f0, ν, 400) 
