using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using ProgressMeter
using OrdinaryDiffEq:Vern9
using ColorSchemes, Colors
include(srcdir("print_fig.jl"))

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
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    ds = CoupledODEs(coupled_lorenz!, zeros(6), [α, β, ε, γ]; diffeq)

    yg = range(-50, 50; length = 400)
    grid = ntuple(x -> yg, 6)
    mapper = AttractorsViaRecurrences(ds, grid; Δt = 0.1, maximum_iterations = Int(1e9),
    consecutive_recurrences = 10000,
    attractor_locate_steps = 1000)
    y1r = range(10, 24, length = res)
    y2r = range(10, 24, length = res)
    bsn = @showprogress [ mapper([1, 1, y1, 1, 1, y2]) for y1 in y1r, y2 in y2r] 
    grid = (y1r,y2r); att = mapper.bsn_nfo.attractors
    return @strdict(bsn, grid, att, res)
end

let res = 1200
α = 10; β = 24.76; ε = 1.1; γ = 8/3; 
params = @strdict res α β ε γ
cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
print_fig(params, "coupled_lorenz", compute_lorenz; force = false, xlab = L"z_1", ylab = L"z_2", cmap) 
# att = get_att(params, "coupled_lorenz", compute_lorenz; force = false) 
# @show att 
end
