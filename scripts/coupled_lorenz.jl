using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using DynamicalSystems
using Attractors
using OrdinaryDiffEq:Vern9

# PHYSICAL REVIEW E 96, 062203 (2017)
# Coupled Lorenz oscillators near the Hopf boundary: Multistability, intermingled basins, and quasiriddling
# Thierry T. Wontchui, Joseph Y. Effa, H. P. Ekobena Fouda, Sangeeta R. Ujjwal, and Ram Ramaswamy
# https://doi.org/10.1103/PhysRevE.96.062203
function coupled_lorenz!(du, u, p, t)
    α, β, ε, γ = p
    x1, y1, z1, x2, y2, z2 = u
    du[1] = α*(y1-x1)
    du[2] = β*x1 - y1 - x1*z1 
    du[3] = -γ*z1 + x1*y1 + ε*(z2 - z1)
    du[4] = α*(y2 - x2)
    du[5] = β*x2 - y2 - x2*z2
    du[6] = -γ*z2 + x2*y2 + ε*(z1 - z2)
end


function compute_lorenz(di::Dict)
    @unpack res, α, β, ε, γ = di
    ds = ContinuousDynamicalSystem(coupled_lorenz!, zeros(6), [α, β, ε, γ], (J,z0, p, n) -> nothing)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    yg = range(-37, 37; length = 10001)
    grid = ntuple(x -> yg, 6)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, Δt = .1,   
        mx_chk_fnd_att = 100, show_progress = true,
        mx_chk_loc_att = 100, mx_chk_safety = Int(1e8), diffeq)

    y1r = range(0, 20, length = res)
    y2r = range(0, 20, length = res)
    bsn = [ mapper([1, 1, y1, 1, 1, y2]) for y1 in y1r, y2 in y2r] 
    grid = (y1r,y2r); att = mapper.bsn_nfo.attractors
    return @strdict(bsn, grid, att, res)
end


function print_fig(w, h, α, β, ε, γ, res)
    params = @strdict res α β ε γ
    data, file = produce_or_load(
        datadir("basins"), params, compute_lorenz;
        prefix = "coupled_lorenz", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    xg,yg = grid
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"\rho", xlabel = L"z",            
            yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, yg, xg, bsn'; rasterize = 1)
    save(string("../plots/coupled_lorenz_", res,".png"),fig)
end


α = 10; β = 24.76; ε = 1.1; γ = 8/3; 
print_fig(600,600, α, β, ε, γ, 1000)
