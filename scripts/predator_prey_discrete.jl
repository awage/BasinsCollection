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
    return nothing
end



function compute_pred_prey(di)
    @unpack m, k, res = di
    df = DeterministicIteratedMap(predator_prey!, rand(2), [m,k]) 
    y1 = collect(range(-1, 20, length = 40000))
    # y1 = range(0, 10, length = 600001)
    pow = 2; x1 = range(0, 20.0^(1/pow); length = 40000).^pow
    grid_rec = (x1, y1)
    mapper = AttractorsViaRecurrences(df, grid_rec,
        consecutive_recurrences = 10000,
        consecutive_attractor_steps = 100
        ,attractor_locate_steps = 10000)
    x = range(0.1, 2, length = res)
    y = range(0.01, 0.2, length = res)
    grid = (x,y)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, res)
end


let res = 1200
m = 0.104; k = 7.935; 
params = @strdict k m res
print_fig(params, "pred_prey", compute_pred_prey; ylab = L"y_0", xlab = L"x_0", force = true)
end
