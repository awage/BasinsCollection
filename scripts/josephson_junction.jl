using DrWatson
@quickactivate "AdvancedBasins" # exports Attractors, GLMakie and other goodies in `src`
using Attractors
using OrdinaryDiffEq
using LaTeXStrings
using CairoMakie
# Cell-mapping techniques applied to the rf-driven Josephson junction
# Article in Physical review A, Atomic, molecular, and optical physics · October 1987
# DOI: 10.1103/PhysRevA.36.2455 · Source: PubMed
function josephson_junction!(du, u, p, t)
    β = p[1]; idc = p[2]; irf = p[3]; Ω = p[4]
    du[1] = u[2]
    du[2] = -β^(-0.5)*u[2] - sin(u[1]) + idc + irf*sin(Ω*t) 
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



function compute_basins_josephson(di::Dict)
    @unpack β, idc, irf, Ω, res = di
    ds = ContinuousDynamicalSystem(josephson_junction!, zeros(2), [β, idc, irf, Ω], (J,z0, p, n) -> nothing)
    condition(u,t,integrator) = (integrator.u[1] < -π  || integrator.u[1] > π)
    cb = DiscreteCallback(condition,affect!)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8, callback = cb)
    pstrob = stroboscopicmap(ds, 2π/Ω; diffeq)
    y1 = range(-pi, pi, length = 8000)
    y2 = range(-5, 5, length = 8000)
    mapper = AttractorsViaRecurrences(pstrob, (y1,y2); sparse = true,    
        mx_chk_fnd_att = 1000,
        mx_chk_loc_att = 1000, safety_counter_max = Int(1e7) )
    y1 = range(-pi, pi, length = res)
    y2 = range(-3, 1, length = res)
    bsn, att = basins_of_attraction(mapper, (y1,y2))
    grid = (y1,y2)
    return @strdict(bsn, att, grid, β, idc, irf, Ω, res)
end


function print_fig(w, h, cmap, β, idc, irf, Ω, res)
    params = @strdict β idc irf Ω res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_josephson;
        prefix = "josephson", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    xg, yg = grid
    fig = Figure(size = (w, h))
    ax = Axis(fig[1,1], ylabel = L"$\dot{x}$", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    if isnothing(cmap)
        heatmap!(ax, xg, yg, bsn, rasterize = 1)
    else
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    end
    save(string(projectdir(), "/plots/basins_josephson",res,".png"),fig)
end

β = 25; idc = 1.878; irf = 10.198; Ω = 1;
print_fig(600,600, nothing, β,  idc, irf,  Ω, 700) 
