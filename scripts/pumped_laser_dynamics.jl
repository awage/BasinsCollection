using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEqVerner
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
    m,f = p
    # m = 0.8; f = 70.2e3
    a = 6.6207e7; b =7.4151e6; c = 0.0163;
    d = 4.0763e3; ρ = 0.3075; 
    pm = 506*(1 + m*sin(2*π*f*t))
    x,y = u
    du[1] = a*x*y - b*x + c*(y + ρ)
    du[2] = -d*x*y -(y +ρ) + pm*(1-exp(-18*(1-(y+ρ)/0.615))) 
    return nothing
end




function compute_basins_lasers(di::Dict)
    @unpack  res, m ,f = di
    ds = CoupledODEs(pumped_laser!, rand(2), [m,f])
    smap = StroboscopicMap(ds, 1/f)
    pow = 2; xg = range(0, 1^(1/pow); length = 5000).^pow
    yg =  collect(range(0, 1; length = 5000))
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(smap, grid;
    consecutive_recurrences = 2000) 
    # consecutive_attractor_steps = 10)
    yr = range(0.09, 0.13, length = res)
    xr = range(0.0, 15, length = res)
    bsn = @showprogress [ mapper([x,y]) for x in xr, y in yr]
    grid = (xr,yr)
    att = extract_attractors(mapper)
    return @strdict(bsn, grid, res, att)
end


let res = 1200
m = 0.8; f = 70.2e3
# m = 1; f = 80e3
params = @strdict res m f
print_fig(params, "basins_lasers", compute_basins_lasers; force = false)
end

