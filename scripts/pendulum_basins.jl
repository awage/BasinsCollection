using DrWatson
@quickactivate
using OrdinaryDiffEq
using DynamicalSystems
using CairoMakie
using LaTeXStrings
using ColorSchemes


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
    df = ContinuousDynamicalSystem(forced_pendulum,rand(2), [d, F, ω])
    condition(u,t,integrator) = (integrator.u[1] < -π  || integrator.u[1] > π)
    cb = DiscreteCallback(condition,affect!)
    diffeq = (reltol = 1e-9,  alg = Vern9(), callback = cb)
    xg = range(-pi,pi,length = res)
    yg = range(-4.,4.,length = res)
    smap = stroboscopicmap(df, 2*pi/ω; diffeq)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    bsn, att = basins_of_attraction(mapper)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, d, F, ω, res)
end


function print_fig(w, h, cmap, d, F, ω, res) 
    # d = 0.2; F = 1.3636363636363635; ω = 0.5 # Parameters for Riddled Basins
    # res = 2000
    params = @strdict d F ω res

    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_pend;
        prefix = "pendulum", storepatch = false, suffix = "jld2", force = false
    )


    @unpack bsn, grid = data
    xg ,yg = grid
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"\dot{\theta}", xlabel = L"\theta", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    save(string("../plots/basins_pendulum_", ω, ".png"),fig)
end

function get_Sb(d, F, ω, res)
    params = @strdict d F ω res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_pend;
        prefix = "pendulum", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    return basin_entropy(bsn)
end
