using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors
using StaticArrays


# Multistability, phase diagrams, and intransitivity in the Lorenz-84 low-order atmospheric circulation model
# Chaos 18, 033121 (2008); https://doi.org/10.1063/1.2953589
@inline @inbounds function lorenz84(u, p, t)
    F = p[1]; G = p[2]; a = p[3]; b = p[4];
	x = u[1]; y = u[2]; z = u[3];
    dx = -y^2 -z^2 -a*x + a*F
    dy = x*y - y - b*x*z +G
	dz = b*x*y + x*z -z
    return SVector{3}(dx, dy, dz)
end


function compute_lorenz84(di)
    @unpack F, G, a, b, res = di
#F=6.846; G=1.287; a=0.25; b=4.;
    ds = Systems.lorenz84(F = F, G = G, a = a, b = b)
    xg=range(-2.,2.,length = 2000)
    yg=range(-2.,2.,length = 2000)
    zg=range(-2.,2.,length = 2000)
    grid_rec = (xg,yg,zg)
    mapper = AttractorsViaRecurrences(ds, grid_rec,
            # mx_chk_lost = 10, 
            mx_chk_fnd_att = 300, 
            mx_chk_loc_att = 300, 
            # mx_chk_att = 2,
             sparse = true
    )
    xg=range(-2.,2.,length = res)
    yg=range(-2.,2.,length = res)
    bsn = [ mapper([x,y,0.]) for x in xg, y in yg]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(bsn, att, grid, F, G, a, b, res)
end



function print_fig(w, h, cmap, F, G, a, b, res) 
    params = @strdict F G a b res

    data, file = produce_or_load(
        datadir("basins"), params, compute_lorenz84;
        prefix = "lorenz84", storepatch = false, suffix = "jld2", force = false
    )


    @unpack bsn, grid = data
    xg ,yg = grid
    fig = Figure(size = (w, h))
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
    save(string("../plots/lorenz84_", res, ".png"),fig)
end



 

F=6.846; G=1.287; a=0.25; b=4.; #res = 600
params = @strdict F G a b res
print_fig(params, "lorenz84", compute_lorenz84; ylab = L"\dot{\theta}", xlab = L"\theta")
