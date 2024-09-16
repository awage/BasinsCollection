using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie

function henon_map!(dz, z, p, n)
    xn = z[1]; yn = z[2]
    μ, j = p
    dz[1] = 1 -μ*xn^2 +yn
    dz[2] = -j*xn  
    return
end


function compute_henon(di::Dict)
    @unpack μ, j, res = di
    ds = DeterministicIteratedMap(henon_map!, [1.0, 0.0], [μ, j])
    yg = xg = range(-2., 2., length = 25000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg))
    yg = xg = range(-2, 2, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, μ, j, res)
end


μ = 1.08
j = 0.9
# res = 1000
params = @strdict res μ j
print_fig(params, "henon", compute_henon)
