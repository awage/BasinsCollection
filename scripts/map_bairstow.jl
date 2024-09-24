using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie

# On noninvertible mappings of the plane: Eruptions
# Lora Billings;
# James H. Curry
# Chaos 6, 108–120 (1996)
# https://doi.org/10.1063/1.166158
function bairstow_map!(dz, z, p, n)
    u = z[1]; v = z[2]
    a = p
    dz[1] = (u^3 + u*(v-a+1) + a)/(2*u^2 + v)
    dz[2] = (v*(u^2 + a - 1) + 2*a*u)/(2*u^2 + v)
    return
end


function compute_bairstow(di::Dict)
    @unpack a, res = di
    ds = DeterministicIteratedMap(bairstow_map!, [1.0, 0.0], a)
    yg = xg = range(-4., 4., length = 25000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg))
    yg = xg = range(-4, 4, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, μ, j, res)
end


# res = 1000
a = 0.8
params = @strdict res a 
print_fig(params, "bairstow_map", compute_bairstow; xlab = L"u", ylab = L"v")
