using DrWatson
@quickactivate
using OrdinaryDiffEq:Vern9
using Attractors
using CairoMakie
using LaTeXStrings

# Eur. Phys. J. Special Topics 226, 1979–1985 (2017)
# https://doi.org/10.1140/epjst/e2017-70037-1
# Megastability: Coexistence of
# a countable infinity of nested attractors
# in a periodically-forced oscillator with
# spatially-periodic damping
# Julien C. Sprott , Sajad Jafari , Abdul Jalil M. Khalaf , and Tomasz Kapitaniak
function megastable(u, p, t)
    x,y = u
    du1 = y
    du2 = -0.33^2*x + y*cos(x) + sin(0.73*t)
    return SVector{2}(du1, du2)
end

function compute_megastable(di::Dict)
    @unpack  res = di
    diffeq = (reltol = 1e-6,  alg = Vern9())
    df = CoupledODEs(megastable,rand(2), []; diffeq)
    xg = yg = range(-5000,5000,length = 50000)
    smap = StroboscopicMap(df, 2π/0.73)
    mapper = AttractorsViaRecurrences(smap, (xg, yg);
    consecutive_recurrences = 1000,
    attractor_locate_steps = 1000)
    xg = yg = range(-20,20,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid,  res)
end


res = 1200
params = @strdict res
print_fig(params, "megastable", compute_megastable; force = false)

