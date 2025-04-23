using DrWatson
@quickactivate
using OrdinaryDiffEqVerner
using Attractors
using CairoMakie
using LaTeXStrings
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))
# Pramana – J. Phys. (2018) 91:11
# https://doi.org/10.1007/s12043-018-1581-6
# Carpet oscillator: A new megastable nonlinear oscillator
# with infinite islands of self-excited and hidden attractors
# YANXIA TANG, HAMID REZA ABDOLMOHAMMADI
# YE TIAN and TOMASZ KAPITANIAK5
# , ABDUL JALIL M KHALAF
# Equations of motion:
function carpet(u, p, t)
    x,y = u
    du1 = sin(0.1*y)
    du2 = -sin(0.1*x) + sin(0.1*y)*cos(x) 
    return SVector{2}(du1, du2)
end

function compute_carpet(di::Dict)
    @unpack  res = di
    diffeq = (reltol = 1e-6,  alg = Vern9())
    df = CoupledODEs(carpet,rand(2), []; diffeq)
    xg = yg = range(-5000,5000,length = 50000)
    mapper = AttractorsViaRecurrences(df, (xg, yg);
    consecutive_recurrences = 1000,
    attractor_locate_steps = 1000)
    xg = yg = range(-100,100,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid,  res)
end


res = 1200
params = @strdict  res
print_fig(params, "carpet", compute_carpet; force = false)


