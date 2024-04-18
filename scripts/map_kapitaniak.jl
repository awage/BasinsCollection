using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9

# Bifurcations from locally to globally riddled basins
# T. Kapitaniak, Yu. Maistrenko, A. Stefanski, and J. Brindley
# https://doi.org/10.1103/PhysRevE.57.R6253
function kapitaniak_map!(dz, z, p, n)
    xn = z[1]; yn = z[2]
    l,pp,d1,d2 = p
    dz[1] = pp*xn + l/2*(1-pp/l)*(abs(xn+1/l)-abs(xn-1/l))+d1*(yn-xn)
    dz[2] = pp*yn + l/2*(1-pp/l)*(abs(yn+1/l)-abs(yn-1/l))+d2*(xn-yn)
    return
end

# dummy function to keep the initializator happy
function kapitaniak_map_J(J, z0, p, n)
    return
end



function compute_kapitaniak(di::Dict)
    @unpack l, pp, d1, d2, res = di
    ds = DeterministicIteratedMap(kapitaniak_map!, [1.0, 0.0], [l, pp, d1, d2])
    yg = xg = range(-3., 3., length = 10000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        consecutive_recurrences = 10000,
        attractor_locate_steps = 10000, maximum_iterations = Int(1e8), show_progress = true)
    yg = xg = range(-2, 2, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end


# l = 1.3; pp = -2.; d1 = 0.725; d2 = 0.725;
l = √2; pp = -√2; d1 = d2 = -0.935;


params = @strdict res l pp d1 d2
print_fig(params, "kapitaniak", compute_kapitaniak) 
