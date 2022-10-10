using DrWatson
@quickactivate "AdvancedBasins" # exports DynamicalSystems, GLMakie and other goodies in `src`
using CairoMakie
using LaTeXStrings


include("kur_halekotte.jl")


function makesim(di::Dict)
    @unpack ni, res = di
    bsn, att, grid = get_basins(ni, res)
    return @strdict(bsn, att, grid, ni, res)
end

function print_fig(ni, res)

    params = @strdict ni res

    data, file = produce_or_load(
        datadir("basins"), params, makesim;
        prefix = "basins_kur", storepatch = false, suffix = "jld2", force = false
    )


    @unpack bsn, grid = data

    xg, yg = grid

    fig = Figure(resolution = (600, 600))
    ax = Axis(fig[1,1], ylabel = L"$\phi$", xlabel = L"\omega", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    # heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    heatmap!(ax, xg, yg, bsn, rasterize = 1)
    save(string("../plots/basins_", ni,".png"),fig)

end


# print_fig(58,10)
