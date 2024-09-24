using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter

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

    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e9)
    ds = CoupledODEs(memristor!, rand(4), [α, β, γ, g, r, m]; diffeq)
    xg = yg = zg = tg = range(-3, 3; length = 10001)
    grid = (xg, yg, zg, tg)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, Δt = 0.01,   
        consecutive_basin_steps = 200,
        consecutive_recurrences = 1000,
        attractor_locate_steps = 1000, maximum_iterations = Int(1e8), show_progress = true)
    yr = range(-2, 2, length = res)
    xr = range(-2, 2, length = res)
    bsn = @showprogress [ mapper([x,y,0.]) for x in xr, y in yr]
    grid = (xr,yr)
    return @strdict(bsn, grid, res)
end

α = 1; β = 0.05; γ = 0.5; g = 0.03; r = 5.8; m = 11; 
params = @strdict α β γ g r m res
print_fig(params, "hidden_memristor", compute_basins_memristor; force = false) 
