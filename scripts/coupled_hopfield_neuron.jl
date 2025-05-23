using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using ProgressMeter
using OrdinaryDiffEqVerner
using ColorSchemes, Colors
include(srcdir("print_fig.jl"))


function coupled_hopfield!(dz, z, p, t)
    x, u, y, z, v, w = z 
    α, β, δ, b1, b2 = p
    f(x) = tanh(2*x - 3) + tanh(2*x + 3) - 2*tanh(2*x)

    # Equations
    dz[1] = u                                             
    dz[2] = -δ * u - (x + α * y) + b1 * (v^2 - 0.5) * f(x + α * y) 
    dz[3] = z                                             
    dz[4] = -δ * z - (y + β * x) + b2 * (w^2 - 0.5) * f(y + β * x) 
    dz[5] = f(x + α * y) - 2 * v                  
    dz[6] = f(y + β * x) - 2 * w                 
end


function compute_hpf_cpld(di::Dict)
    @unpack res, α, β, δ, b1, b2 = di
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    ds = CoupledODEs(coupled_hopfield!, zeros(6), [α, β, δ, b1, b2]; diffeq)

    yg = range(-10, 10; length = 1001)
    grid = ntuple(x -> yg, 6)
    mapper = AttractorsViaRecurrences(ds, grid; Δt = 1., maximum_iterations = Int(1e9),
    consecutive_recurrences = 1000,
    # Ttr = 1000,
    attractor_locate_steps = 1000,
    consecutive_attractor_steps = 4)
    y1r = range(-2, 2, length = res)
    y2r = range(-2, 2, length = res)
    bsn = @showprogress [ mapper([y1, 1., y2, 1., 1., 5.]) for y1 in y1r, y2 in y2r] 
    grid = (y1r,y2r); att = mapper.bsn_nfo.attractors
    return @strdict(bsn, grid, att, res)
end

res = 300
δ = 1.5; α = -0.1; b1 = 3.0; b2 = 3.0; β = -α
params = @strdict res δ α b1 b2 β
# cmap = ColorScheme([ RGB(0.9,0.2,0.1), RGB(1,1,1), RGB(0,0,0) ] )
print_fig(params, "hopfield_neuron", compute_hpf_cpld; force = true, xlab = L"z_1", ylab = L"z_2") 
