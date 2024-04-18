using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors
using ChaosTools
using StaticArrays

# https://doi.org/10.48550/arXiv.2211.06921
# Disipative nontwist system
function dsnm!(dz, z, p, n)
    x,y = z
    γ = 0.1; 
    a, b = p
    dz[2] = (1-γ)*y - b*sin(2π*x)
    dz[1] = rem(x + a*(1-dz[2]^2),1, RoundNearest)
end



function compute_dsnm(di)
    @unpack a, b, res = di
    u0 = [0.; 0.]
    p = [a, b]
    df = DiscreteDynamicalSystem(dsnm!, u0, p) 
    x1 = range(-2, 2, length = 10001)
    y1 = range(-5, 5, length = 10001)
    grid_rec = (x1, y1)
    mapper = AttractorsViaRecurrences(df, grid_rec,
            # mx_chk_lost = 10, 
            mx_chk_fnd_att = 1000, 
            mx_chk_loc_att = 1000, 
            mx_chk_att = 4,
            safety_counter_max = Int(1e8),
            sparse = true, Ttr = 300)
    x = range(0., 1, length = res)
    y = range(-2, 2, length = res)
    grid = (x,y)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, res)
end
 
a = 0.55; b = 0.45; res = 1000
params = @strdict a b res
print_fig(params, "dsnm", compute_dsnm)

