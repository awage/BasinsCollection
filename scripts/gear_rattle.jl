using DrWatson
@quickactivate
using OrdinaryDiffEqVerner
using Attractors
using CairoMakie
using LaTeXStrings
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))

# International Journal of Bifurcation and Chaos, Vol. 19, No. 1 (2009) 203–224
# BASINS OF ATTRACTION IN NONSMOOTH MODELS OF GEAR RATTLE
# JOANNA F. MASON
# PETRI T. PIIROINEN
# R. EDDIE WILSON and MARTIN E. HOMER
function B(x,β) 
    if x ≥ β 
        return(x - β)
    elseif abs(x) < β
        return(0)
    elseif x ≤ -β
        return(x + β)
    end
end

# Equations of motion:
function gear_rattle(u, p, t)
    @inbounds begin
    β, δ, ε = p
    κ = 100
    du1 = u[2]
    du2 = -δ*u[2] -2*κ*B(u[1],β) + 4*π*δ -4*π^2*ε*cos(2π*t) - 2*π*ε*δ*sin(2π*t)
    return SVector{2}(du1, du2)
    end
end


function compute_gear_rattle(di::Dict)
    @unpack δ, ε, β ,res = di

    diffeq = (reltol = 1e-10,  alg = Vern9(),)
    df = CoupledODEs(gear_rattle, rand(2), [ β, δ, ε]; diffeq)
    xg = range(-4.,4.,length = 1001)
    yg = range(-4.,4.,length = 1001)
    smap = StroboscopicMap(df, 2*pi)
    mapper = AttractorsViaRecurrences(smap, (xg, yg);
            consecutive_attractor_steps = 2,
            # consecutive_basin_steps = 400,
            consecutive_recurrences = 400,
            attractor_locate_steps = 200)
    xg = range(-β+1e-4,β-1e-4,length = res)
    yg = range(-4.,4.,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


β = 0.6; δ = 0.6; ε = 0.1; res = 1200
params = @strdict β δ ε res
print_fig(params, "gear_rattle", compute_gear_rattle; ylab= L"\dot{\theta}", xlab= L"\theta", force = false)
# att = get_att(params, "gear_rattle", compute_gear_rattle)

# function gear_rattle!(du,u, p, t)
#     δ, ε = p
#     κ = 100
#     du[1] = u[2]
#     du[2] = -δ*u[2] + 4*π*δ -4*π^2*ε*cos(2π*t) - 2*π*ε*δ*sin(2π*t)
# end
