using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
# From attractor to chaotic saddle: a tale of transverse  instability
# Peter Ashwin, Jorge Buescu and Ian Stewart
# Nonlinearity 9 (1996) 703–737. 
# https://doi.org/10.1088/0951-7715/9/3/006
function cplog_ashwin(dz, z, p, n)
    x1 = z[1]; x2 = z[2]
    α, ν, ε = p
    dz[1] = 3*√3/2*x1*(x1^2-1) + ε*x1*x2^2
    dz[2] = ν*exp(-α*x1^2)*x2 + x2^3
    return
end


function compute_ashwin(di::Dict)
    @unpack α, ν, ε, res = di
    ds = DeterministicIteratedMap(cplog_ashwin, [1.0, 0.0], [α, ν, ε])
    yg = xg = range(-20., 20., length = 25000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        consecutive_recurrences = 1000,
        attractor_locate_steps = 1000, maximum_iterations = Int(1e8), show_progress = true)
    yg = range(-0.7, 0.7, length = res)
    xg = range(-1.5, 1.5, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid,  res)
end



# a = 3.57480493875920; ε = -0.2
ν = 1.28; α = 0.7; ε = 0.5; 
params = @strdict res α ε ν
print_fig(params, "cpashwin", compute_ashwin) 

