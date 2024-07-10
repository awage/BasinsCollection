using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter
include(srcdir("print_fig.jl"))
# using Plots


# Control of basins of attraction in a multistable fiber laser
# A.N. Pisarchik a,∗ , R. Jaimes-Reategui b
# Physics Letters A 374 (2009) 228–234
# doi:10.1016/j.physleta.2009.10.061

# Model has been normalized see https://doi.org/10.3390/photonics11020176
# for more details

function pumped_laser!(du, u, p, t)
    a = 6.6207e7; b =7.4151e6; c = 0.0163;
    d = 4.0763e3; ρ = 0.3075; m = 0.8; f = 70.2e3
    pm = 506*(1 + m*sin(2*π*f*t))
    x,y = u
    du[1] = a*x*y - b*x + c*(y + ρ)
    du[2] = -d*x*y -(y +ρ) + pm*(1-exp(-18*(1-(y+ρ)/0.615))) 
    return nothing
end




function compute_basins_lasers(di::Dict)
    @unpack  res = di
    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e9)
    ds = CoupledODEs(pumped_laser!, rand(2), []; diffeq)
    pow = 2; xg = range(0, 50^(1/pow); length = 20000).^pow
    # yg =  collect(range(0, 1; length = 20000))
    yg =  collect(range(0, 50; length = 20000))
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid; Δt = 1e-7, 
    consecutive_recurrences = 1000, 
    consecutive_attractor_steps = 20)
    yr = range(0.09, 0.13, length = res)
    xr = range(0.0, 15, length = res)
    bsn = @showprogress [ mapper([x,y]) for x in xr, y in yr]
    grid = (xr,yr)
    att = extract_attractors(mapper)
    return @strdict(bsn, grid, res, att)
end


let res = 1200
params = @strdict res
print_fig(params, "basins_lasers", compute_basins_lasers; force = false)
end

# diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e9)
# ds = CoupledODEs(pumped_laser!, rand(2), []; diffeq)
# y,t = trajectory(ds, 0.006, [0.01*rand(),rand()];Δt  = 0.000001); 

# lines(y[1:1000,1])

