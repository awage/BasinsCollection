using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using DynamicalSystems
using OrdinaryDiffEq:Vern9
using ProgressMeter
# using Plots

# Multistable dynamics and control of a new 4D memristive chaotic
# Sprott B system
# Ramesh Ramamoorthy a, Karthikeyan Rajagopal b, Gervais Dolvis Leutcho c,d,
# Ondrej Krejcar e,f, Hamidreza Namazi e,g,∗, Iqtadar Hussain h
# https://doi.org/10.1016/j.chaos.2022.111834
function memristor!(du, u, p, t)
    α, β, γ, g, r, m = p
    x,y,z,t = u
    a = 1 
    W(u) = α + γ*abs(u) + β*u^2
    du[1] = r*y*z + g
    du[2] = x - y
    du[3] = 1 - m*W(t)*x*y
    du[4] = a*x*y - t
end




function compute_basins_memristor(di::Dict)
    @unpack α, β, γ, g, r, m, res = di

    ds = ContinuousDynamicalSystem(memristor!, rand(4), [α, β, γ, g, r, m])
    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e9)
    xg = yg = zg = tg = range(-3, 3; length = 10001)
    grid = (xg, yg, zg, tg)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, Δt = 0.01,   
        mx_chk_hit_bas = 200,
        mx_chk_fnd_att = 1000,
        mx_chk_loc_att = 1000, safety_counter_max = Int(1e8), diffeq)
    yr = range(-2, 2, length = res)
    xr = range(-2, 2, length = res)
    bas = @showprogress [ mapper([x,y,0.]) for x in xr, y in yr]
    grid = (xr,yr)
    return @strdict(bas, grid, res)
end


function print_fig(w,h,cmap, α, β, γ, g, r, m, res)
    params = @strdict α β γ g r m res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_memristor;
        prefix = "hidden_memristor", storepatch = false, suffix = "jld2", force = true
    )
    @unpack bas, grid = data
    xg, yg = grid

    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"y", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    if isnothing(cmap)
        heatmap!(ax, xg, yg, bas, rasterize = 1)
    else
        heatmap!(ax, xg, yg, bas, rasterize = 1, colormap = cmap)
    end

    save(string("../plots/basins_hidden_memristor", res,".png"),fig)

end

α = 1; β = 0.05; γ = 0.5; g = 0.03; r = 5.8; m = 11; 
print_fig(600, 500, nothing, α, β, γ, g, r, m, 500) 

