using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes

function B(x) 
    β = 0.6
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
    # d = p[1]; F = p[2]; omega = p[3]
    δ, ε = p
    κ = 100
    du1 = u[2]
    du2 = -δ*u[2] -2*κ*B(u[1]) + 4*π*δ -4*π^2*ε*cos(2π*t) - 2*π*ε*δ*sin(2π*t)
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


function compute_gear_rattle(di::Dict)
    @unpack δ, ε, res = di
    # condition(u,t,integrator) = (integrator.u[1] < -β  || integrator.u[1] > β)
    # cb = DiscreteCallback(condition,affect!)
    # diffeq = (reltol = 1e-9,  alg = Vern9(), callback = cb)
    diffeq = (reltol = 1e-8,  alg = Vern9(),)
    df = CoupledODEs(gear_rattle, rand(2), [δ, ε]; diffeq)
    xg = range(-4.,4.,length = 10001)
    yg = range(-4.,4.,length = 10001)
    mapper = AttractorsViaRecurrences(df, (xg, yg);
            consecutive_basin_steps = 200,
            consecutive_recurrences = 1000,
            attractor_locate_steps = 2000)
    xg = range(-1.,1.,length = res)
    yg = range(-4.,4.,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


δ = 0.6; ε = 0.1; res = 300
params = @strdict δ ε res
print_fig(params, "gear_rattle", compute_gear_rattle; ylab= L"\dot{\theta}", xlab= L"\theta", force = true)
# att = get_att(params, "gear_rattle", compute_gear_rattle)
