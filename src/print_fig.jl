# using JLD2
using LaTeXStrings
using CairoMakie

function print_fig(params, sys_name, fun_name; w = 800, h = 800, cmap = nothing, xlab = L"x", ylab = L"y", force = false, format = :pdf)
    println(sys_name)
    data, file = produce_or_load(
        datadir("basins"), params, fun_name;
        prefix = sys_name, storepatch = false, suffix = "jld2", force = force
    )
    @unpack bsn, grid = data
    xg, yg = grid


    if format == :pdf
        fig = Figure(figure_padding = 23, size = (w, h))
        ax = Axis(fig[1,1], ylabel = ylab, xlabel = xlab, 
                yticklabelsize = 30, 
                xticklabelsize = 30, 
                ylabelsize = 40, 
                xlabelsize = 40, 
                xticklabelfont = "NewComputerModern", 
                yticklabelfont = "NewComputerModern")
    else
        fig = Figure(size = (200, 200))
        ax = Axis(fig[1,1], 
            xticksvisible = false, 
            yticksvisible = false, 
            xticklabelsvisible = false, 
            yticklabelsvisible = false)
    end

    if isnothing(cmap)
        println("Number of basins: " , length(unique(bsn)))
        cmap = :mk_12
        if -1 ∈ unique(bsn)
            bsn = Float32.(bsn)
            bsn[bsn .== -1] .= NaN
        end
        # cmap = :berlin10
        # cmap = :seaborn_bright
        # cmap = :dracula
        # cmap = :flag
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    else
        println("Number of basins: " , length(unique(bsn)))
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    end

    if format == :pdf
        resize_to_layout!(fig)
        save(plotsdir(savename(sys_name,params,"pdf")),fig)
    else
        save(plotsdir(savename(sys_name,params,"png")),fig)
        println("![",sys_name,"](./plots/", savename(sys_name,params,"png"),")")
    end
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

