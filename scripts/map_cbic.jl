using DrWatson
@quickactivate
using Attractors
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
include("../src/print_fig.jl")

# Wada bifurcations and partially Wada basin boundaries in a two-dimensional cubic map Yongxiang Zhang , Guanwei Luo 
# http://dx.doi.org/10.1016/j.physleta.2013.03.027
function cubic_map!(dz, z, p, n)
    xn = z[1]; yn = z[2]
    μ, j = p
    dz[1] = yn 
    dz[2] =μ*yn - yn^3 - j*xn  
    return
end

# dummy function to keep the initializator happy
function cubic_map_J(J, z0, p, n)
    return
end


function compute_cbic_map(di::Dict)
    @unpack μ, j, res = di
    ds = DeterministicIteratedMap(cubic_map!, [1.0, 0.0], [μ, j])
    yg = xg = range(-2., 2., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        consecutive_recurrences = 1000,
        attractor_locate_steps = 1000, maximum_iterations = Int(1e8), show_progress = true)
    yg = xg = range(-2, 2, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, μ, j, res)
end

μ = 2.9
j = 0.66
# res = 1000
params = @strdict res μ j
print_fig(params, "cbic_map", compute_cbic_map) 
