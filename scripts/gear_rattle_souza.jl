using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes

# De Souza, S. L. T., & Caldas, I. L. (2001). Basins of Attraction and Transient Chaos in a Gear-Rattling Model. Journal of Vibration and Control, 7(6), 849–862. doi:10.1177/107754630100700605 
# Equations of motion:
function gear_rattle2(u, p, t)
    γ = 0.1; α = 0.48; β = 0.1; ω = 1;
    du1 = u[2]
    du2 = -β*u[2] + γ  -α*β*ω*cos(ω*t) + α*ω^2*sin(ω*t)
    return SVector{2}(du1, du2)
end

function affect!(integrator)
    r = 0.9
    uu = integrator.u
    # @show integrator
    if (round(uu[1]) == 0. && uu[2] > 0) || (round(uu[1]) == -1. && uu[2] < 0) 
        # println("MODIF!")
        set_state!(integrator, SVector(uu[1], -r*uu[2]))
        u_modified!(integrator, true)
    end
    # @show integrator
end

    condition2(u,t,integrator) = (u[1] + 1)*u[1]

function compute_gear_rattle2(di::Dict)
    @unpack res = di
    cb2 = ContinuousCallback(condition2,affect!; abstol = 1e-16, save_positions=(true,true), interp_points = 30)
    diffeq = (reltol = 0, abstol= 1e-14, alg = Vern9(), callback = cb2)
    df = CoupledODEs(gear_rattle2, rand(2); diffeq)
    xg = range(-0.999,-0.001,length = 5001)
    yg = range(-4.,4.,length = 5001)
    mapper = AttractorsViaRecurrences(df, (xg, yg);
            # consecutive_attractor_steps = 4,
            # consecutive_basin_steps = 400,
            consecutive_recurrences = 1000, Ttr = 2)
            # attractor_locate_steps = 2000)
    xg = range(-0.99,-0.02,length = res)
    yg = range(-1,1,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end

res = 1200
params = @strdict res
print_fig(params, "gear_rattle2", compute_gear_rattle2; ylab= L"\dot{s}", xlab= L"s", force = false)
# att = get_att(params, "gear_rattle2", compute_gear_rattle2)

# cb2 = ContinuousCallback(condition2,affect!; abstol = 1e-16, save_positions=(true,true))
# diffeq = (reltol = 0, abstol= 1e-16, alg = Vern9(), callback = cb2)
# df = CoupledODEs(gear_rattle2, rand(2); diffeq)
# y,t = trajectory(df, 100, [-0.1,1])


