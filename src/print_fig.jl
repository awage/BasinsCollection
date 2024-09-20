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
    ax = Axis(fig[1,1], ylabel = ylab, xlabel = xlab, 
            yticklabelsize = 40, 
            xticklabelsize = 40, 
            ylabelsize = 40, 
            xlabelsize = 40, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    # ax = Axis(fig[1,1], 
    #     xticksvisible = false, 
    #     yticksvisible = false, 
    #     xticklabelsvisible = false, 
    #     yticklabelsvisible = false)

    if isnothing(cmap)
        println("Number of basins: " , length(unique(bsn)))
        cmap = :mk_12
        # cmap = :berlin10
        # cmap = :seaborn_bright
        # cmap = :dracula
        # cmap = :flag
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    else
        println("Number of basins: " , length(unique(bsn)))
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    end
    save(plotsdir(savename(sys_name,params,"pdf")),fig)
end


function get_att(params, sys_name, fun_name; force = false)
    println(sys_name)
    data, file = produce_or_load(
        datadir("basins"), params, fun_name;
        prefix = sys_name, storepatch = false, suffix = "jld2", force = force
    )
    @unpack bsn, grid, att = data
    return att 
end

