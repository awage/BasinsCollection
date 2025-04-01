using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))


# Equations of motion:
function forced_pendulum(u, p, t)
    @inbounds begin
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] - sin(u[1])+ F*cos(omega*t)
    return SVector{2}(du1, du2)
    end
end

# We have to define a callback to wrap the phase in [-π,π]
function affect!(integrator)
    uu = integrator.u
    if integrator.u[1] < 0
        set_state!(integrator, SVector(uu[1] + 2π, uu[2]))
        u_modified!(integrator, true)
    else
        set_state!(integrator, SVector(uu[1] - 2π, uu[2]))
        u_modified!(integrator, true)
    end
end


function compute_basins_pend(di::Dict)
    @unpack d, F, ω, res = di
    condition(u,t,integrator) = (integrator.u[1] < -π  || integrator.u[1] > π)
    cb = DiscreteCallback(condition,affect!)
    diffeq = (reltol = 1e-9,  alg = Vern9(), callback = cb)
    df = CoupledODEs(forced_pendulum,rand(2), [d, F, ω]; diffeq)
    xg = range(-pi,pi,length = res)
    yg = range(-4.,4.,length = res)
    smap = StroboscopicMap(df, 2*pi/ω)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, d, F, ω, res)
end


let
d = 0.2; F = 1.3636363636363635; ω = 0.5 # Parameters for Riddled Basins
# res = 1000
cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
params = @strdict d F ω res
# print_fig(600, 500, cmap, d, F, ω, 1000)
print_fig(params, "pendulum", compute_basins_pend; cmap, ylab= L"\dot{\theta}", xlab= L"\theta")

d = 0.2; F = 1.66; ω = 1. # Parameters for Wada Basins
params = @strdict d F ω res
cmap = ColorScheme([RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )
print_fig(params, "pendulum", compute_basins_pend; cmap, ylab= L"\dot{\theta}", xlab= L"\theta")
end
