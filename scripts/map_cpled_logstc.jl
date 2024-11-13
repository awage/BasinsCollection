using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
include(srcdir("print_fig.jl"))

# Transverse instability and riddled basins in a system of two coupled logistic maps
# Yu. L. Maistrenko, V. L. Maistrenko, and A. Popovich E. Mosekilde
#https://doi.org/10.1103/PhysRevE.57.2713
function cplog_map(dz, z, p, n)
    xn = z[1]; yn = z[2]
    a, ε = p
    dz[1] = a*xn*(1-xn) + ε*(yn-xn) 
    dz[2] =  a*yn*(1-yn) + ε*(xn-yn)   
    return
end


function compute_cplog(di::Dict)
    @unpack a, ε, res = di
    ds = DeterministicIteratedMap(cplog_map, [1.0, 0.0], [a, ε])
    yg = xg = range(-2., 2., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg);     
        consecutive_recurrences = 100, Ttr = 5000)
    yg = xg = range(-0.05, 1.2, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end


res = 1200
# a = 3.57480493875920; ε = -0.2
a = 3.6; ε = -1.
params = @strdict res a ε
print_fig(params, "cplog", compute_cplog; force = false) 
att = get_att(params, "cplog", compute_cplog) 

