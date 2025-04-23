using DrWatson
@quickactivate # exports Attractors, GLMakie and other goodies in `src`
using Attractors
using OrdinaryDiffEqVerner
using CairoMakie
using LaTeXStrings
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))


@inline @inbounds function duffing(u, p, t)
    d = p[1]; F = p[2]; omega = p[3]
    du1 = u[2]
    du2 = -d*u[2] + u[1] - u[1]^3 + F*sin(omega*t)
    return SVector{2}(du1, du2)
end



function compute_basins_duffing(di::Dict)
    @unpack d, F, ω, res = di
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    ds = CoupledODEs(duffing, rand(2), [d, F, ω]; diffeq)
    xg = yg = range(-2.2,2.2,length = res)
    smap = StroboscopicMap(ds, 2*pi/ω)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, d, F, ω, res)
end


cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )

# res = 600; 
# d = 0.1; F=0.1; ω=0.1;  # smooth boundary
# params = @strdict d F ω res
# print_fig(params, string("duffing_",d), compute_basins_duffing; ylab = L"$\dot{x}$", xlab = L"x", cmap) 

# d = 0.4; F=0.1; ω=0.1;  # smooth boundary
# params = @strdict d F ω res
# print_fig(params, string("duffing_",d), compute_basins_duffing; ylab = L"$\dot{x}$", xlab = L"x", cmap) 

let
    d = 0.05; F=0.098; ω=1.15;  # smooth boundary
    params = @strdict d F ω res
    print_fig(params, string("duffing_",d), compute_basins_duffing; ylab = L"$\dot{x}$", xlab = L"x") 
end
