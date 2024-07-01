using LaTeXStrings
using CairoMakie


function print_fig(params, sys_name, fun_name; w = 1200, h = 1200, cmap = nothing, xlab = L"x", ylab = L"y", force = false)
    println(sys_name)
    data, file = produce_or_load(
        datadir("basins"), params, fun_name;
        prefix = sys_name, storepatch = false, suffix = "jld2", force = force
    )
    @unpack bsn, grid = data
    xg, yg = grid
    fig = Figure(size = (w, h))
    # ax = Axis(fig[1,1], ylabel = ylab, xlabel = xlab, 
    #         yticklabelsize = 30, 
    #         xticklabelsize = 30, 
    #         ylabelsize = 30, 
    #         xlabelsize = 30, 
    #         xticklabelfont = "cmr10", 
    #         yticklabelfont = "cmr10")
    ax = Axis(fig[1,1], 
        xticksvisible = false, 
        yticksvisible = false, 
        xticklabelsvisible = false, 
        yticklabelsvisible = false)
    cmap = :mk_12

    if isnothing(cmap)
        heatmap!(ax, xg, yg, bsn, rasterize = 1)
    else
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    end
    save(plotsdir(savename(sys_name,params,"png")),fig)
end


