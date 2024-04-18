using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors
using StaticArrays

# Organized periodic structures and coexistence of triple attractors in a
# predator–prey model with fear and refuge
# Shilpa Garai a ,1 , N.C. Pati b ,1 , Nikhil Pal a ,∗,1 , G.C. Layek b ,1
# Chaos, Solitons and Fractals 165 (2022) 112833
function predator_prey!(dz, z, p, n)
    x,y = z
    R = 3.2; P = 0.1; A = 2 ; B = 5; c = 0.9; D1 = 0.3; D2 = 0.1; 
    m, k = p
    dz[1] = x*exp(R*m + R*(1-m)/(1+k*y) - D1 - P*x - A*(1-m)*y/(B+(1-m)*x))
    dz[2] = y*exp(c*A*(1-m)*x/(B+(1-m)*x) - D2)
end



function compute_pred_prey(di)
    @unpack m, k, res = di
    u0 = [0.; 0.]
    p = [m, k]
    df = CoupledODEs(predator_prey!, u0, p) 
    x1 = range(-1, 20, length = 600001)
    y1 = range(-1, 30, length = 600001)
    grid_rec = (x1, y1)
    mapper = AttractorsViaRecurrences(df, grid_rec,
        consecutive_recurrences = 1000,
        attractor_locate_steps = 1000, maximum_iterations = Int(1e8), show_progress = true)
    x = range(0.1, 2, length = res)
    y = range(0.01, 0.2, length = res)
    grid = (x,y)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, res)
end



function print_fig(w, h, cmap, k, m, res) 
    params = @strdict k m res

    data, file = produce_or_load(
        datadir("basins"), params, compute_pred_prey;
        prefix = "pred_prey", storepatch = false, suffix = "jld2", force = false
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
    save(string("../plots/pred_prey_discrt", res, ".png"),fig)
end


 
m = 0.104; k = 7.935 
print_fig(600,600, nothing, k, m, 600) 

# p = [m, k]
# df = DiscreteDynamicalSystem(predator_prey!, rand(2), p) 
# u = trajectory(df, 1000, [0.3, 0.2])

