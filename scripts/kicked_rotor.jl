using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors

# Reference: Map with more than 100 coexisting low-period periodic attractors, Ulrike Feudel,  Celso Grebogi, Brian R. Hunt, and James A. Yorke
# PHYSICAL REVIEW E, VOLUME 54, NUMBER 1 1996
function kicked_feudel!(dz, z, p, n)
    xn = z[1]; yn = z[2]
    ν = p[1]; f0 = p[2]
    dz[1] = mod2pi(xn + yn)
    dz[2] = (1 - ν)*yn + f0*sin(xn + yn)
    return
end


function compute_kicked_rotor(di)
    @unpack f0, ν, res = di
    u0 = [0., 0.6]
    df = DeterministicIteratedMap(kicked_feudel!, u0, [ν,f0]) 
    xg = range(0, 2π, length = 10000)
    yg = range(-220, 220, length = 10000)
    grid_rec = (xg, yg)
    mapper = AttractorsViaRecurrences(df, grid_rec,
        consecutive_recurrences = 1000,
        attractor_locate_steps = 1000, maximum_iterations = Int(1e8), show_progress = true)
    xg = range(0, 2π, length = res)
    yg = range(-π, π, length = res)
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, ν, f0 , res)
end


res = 500
f0 = 4.
ν = 0.02
params = @strdict f0 ν res
print_fig(params, "kicked_rotor", compute_kicked_rotor; ylab = L"\dot{\theta}", xlab = L"\theta")
# print_fig(600,600, nothing, f0, ν, 400) 
