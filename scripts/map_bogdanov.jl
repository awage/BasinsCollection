using DrWatson
@quickactivate
using Attractors
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEqVerner
using Colors
using ColorSchemes
include(srcdir("print_fig.jl"))

# International Journal of Bifurcation and ChaosVol. 03, No. 04, pp. 803-842 (1993) Tutorials and ReviewsNo Access
# THE BOGDANOV MAP: BIFURCATIONS, MODE LOCKING, AND CHAOS IN A DISSIPATIVE SYSTEM
# DAVID K. ARROWSMITH, JULYAN H. E. CARTWRIGHT, ALEXIS N. LANSBURY, and COLIN M. PLACE
# https://doi.org/10.1142/S021812749300074X
function bogdanov_map!(dz, z, p, n)
    x, y = z
    μ, k, ε = p
    dz[2] = y + ε*y + k*x*(x-1) + μ*x*y
    dz[1] = x + dz[2]
    return
end


function compute_bogdanov(di::Dict)
    @unpack μ, k, ε, res = di
    ds = DeterministicIteratedMap(bogdanov_map!, [1.0, 0.0], [μ, k, ε])
    yg = xg = range(-2., 2., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); sparse = true,    
        mx_chk_fnd_att = 1000,
        mx_chk_loc_att = 1000, maximum_iterations = Int(1e7), show_progress = true)
    yg = range(-0.8, 0.8, length = res)
    xg = range(-0.6, 1., length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end

μ = -0.1
k = 1.2
ε = 0.0125 
res = 2000
params = @strdict res μ k ε
cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(0.7,0.7,0.7), RGB(1,0,0)] )
print_fig(params, "bogdanov", compute_bogdanov; force = false,  cmap) 
