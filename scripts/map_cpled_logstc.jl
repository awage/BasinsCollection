using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie

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

# dummy function to keep the initializator happy
function cplog_map_J(J, z0, p, n)
    return
end



function compute_cplog(di::Dict)
    @unpack a, ε, res = di
    ds = DeterministicIteratedMap(cplog_map, [1.0, 0.0], [a, ε])
    yg = xg = range(-2., 2., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        consecutive_recurrences = 10000,
        attractor_locate_steps = 10000, maximum_iterations = Int(1e8), show_progress = true)
    yg = xg = range(-0.05, 1.3, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, μ, j, res)
end



# a = 3.57480493875920; ε = -0.2
a = 3.6; ε = -1.

params = @strdict res a ε
print_fig(params, "cplog", compute_cplog) 

