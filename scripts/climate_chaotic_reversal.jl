using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using ColorSchemes, Colors
using Attractors
using StaticArrays
using ProgressMeter
using OrdinaryDiffEq:Vern9
include(srcdir("print_fig.jl"))

 # C. Gissinger
 # "A new deterministic model for chaotic reversals"
 # DOI: 10.1140/epjb/e2012-20799-5
@inline @inbounds function climate_reversal(u, p, t)
    μ, ν, Γ = p
    Q,D,V = u
    dQ = μ*Q  - V*D
    dD = -ν*D + V*Q
    dV = Γ - V + Q*D
    return SVector{3}(dQ, dD, dV)
end

function compute_climate_reversal(di)
    @unpack μ, ν, Γ , res = di
    diffeq = (alg = Vern9(), reltol = 1e-7, maxiters = 1e8)
    ds =  CoupledODEs(climate_reversal, rand(3), [μ, ν, Γ]; diffeq)
    xg = yg = zg = range(-4, 4,length = 250)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); 
                consecutive_recurrences = 3000, 
                attractor_locate_steps = 2000,
                consecutive_attractor_steps = 20,
                )
    xg = range(-2,2,length = res)
    yg = range(-2,2,length = res)
    bsn = @showprogress [ mapper([x,y,-1.]) for x in xg, y in yg]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end

let
    μ = 0.1193; ν = 0.1; Γ = 0.9; res = 1200
    params = @strdict  res μ ν Γ 
    cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(0.7,0.7,0.7), RGB(1,0,0)] )
    print_fig(params, "climate_reversal", compute_climate_reversal; cmap, xlab = L"Q", ylab = L"D", force = false)
end
