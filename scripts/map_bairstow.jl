using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
include(srcdir("print_fig.jl"))

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
    yg = xg = range(-10., 10., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg); 
        maximum_iterations = 200)
    yg = xg = range(-4, 4, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end


res = 1200
a = 0.8
# a = 0.15
params = @strdict res a 
print_fig(params, "bairstow_map", compute_bairstow; xlab = L"u", ylab = L"v", force = true)
att = get_att(params, "bairstow_map", compute_bairstow)
