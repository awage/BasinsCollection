using DrWatson
@quickactivate # exports Attractors, GLMakie and other goodies in `src`
using Attractors
using OrdinaryDiffEq:Vern9
using CairoMakie
using LaTeXStrings
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))

# Nonlinear dynamics of a single-gap terahertz
# split-ring resonator under electromagnetic
# radiation
# Cite as: Chaos 33, 103131 (2023); doi: 10.1063/5.0157489
# Gervais Dolvis Leutcho,
# Lyne Woodward,
# and François Blanchard

@inline @inbounds function split_ring_res(u, p, t)
    σ = 0.38
    β = 0.4
    η = 0.08
    # μ = 35
    # ω = 1.0285
    ω = p[1]
    μ = p[2]
    du1 = u[2]
    du2 = -σ*u[2] - u[1]  + β*u[1]^2 - η*u[1]^3 + μ*cos(ω*t)
    return SVector{2}(du1, du2)
end


function compute_split_ring(di::Dict)
    @unpack  res, ω, μ = di
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    ds = CoupledODEs(split_ring_res, rand(2), [ω, μ]; diffeq)
    xg = yg = range(-50,50, length = 10000)
    smap = StroboscopicMap(ds, 2*pi/ω)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    xg = range(-10,10, length = res)
    yg = range(-15,15, length = res)
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


# cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )

res = 1200; ω = 1.0285; μ = 35 
params = @strdict res ω μ
cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
print_fig(params, "split_ring_res", compute_split_ring; force = false, xlab = L"q_0", ylab = L"i_0", cmap) 

# d = 0.4; F=0.1; ω=0.1;  # smooth boundary
# params = @strdict d F ω res
# print_fig(params, string("duffing_",d), compute_basins_duffing; ylab = L"$\dot{x}$", xlab = L"x", cmap) 
