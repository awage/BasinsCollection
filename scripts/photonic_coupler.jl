using DrWatson
@quickactivate 
using Attractors
using OrdinaryDiffEq:Vern9
using CairoMakie
using LaTeXStrings

# Enhanced stability, bistability, and exceptional points in saturable active photonic couplers
# Yertay Zhiyenbayev, Yannis Kominis, Constantinos Valagiannopoulos, Vassilios Kovanis and Anastasios Bountis
# PHYSICAL REVIEW A 100, 043834 (2019)
# https://doi.org/10.1103/PhysRevA.100.043834
@inline @inbounds function phot_coupler(u, p, t)
    α, β, ε, k  = p
    α1 = 1; β1 = 1; γ = 1.
    k = k*β1
    α2 = -α*α1
    β2 = β*β1 
    A1, A2, φ = u 
    du1 = -α1*A1 - k/2*A2*sin(φ)
    du2 = -α2/(1+ε*A2^2)*A2 + k/2*A1*sin(φ)
    du3 = (β2 - β1) + γ*(A2^2 - A1^2) + k/2*(A1/A2 - A2/A1)*cos(φ)
    return SVector{3}(du1, du2, du3)
end


function compute_phot_coupler(di::Dict)
    @unpack  res, α, β, ε, k = di
    diffeq = (;reltol = 1e-6, alg = Vern9(), maxiters = 1e6)
    ds = CoupledODEs(phot_coupler, rand(3), [α, β, ε, k]; diffeq)
    xg = yg = range(0,50, length = 10000)
    φ  = range(-50,50, length = 10000)
    mapper = AttractorsViaRecurrences(ds, (xg, yg, φ);
    attractor_locate_steps = 1000)
    xg = range(0.1, 10, length = res)
    yg = range(0.1, 10, length = res)
    grid = (xg,yg)
    # bsn, att = basins_of_attraction(mapper, grid)
    bsn = @showprogress [ mapper([x,y, π]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid, res)
end


let res = 250; 
α = 2; β = 1.5; ε = 0.5; k = 5 
params = @strdict res α β ε k 
print_fig(params, "photonic_coupler", compute_phot_coupler; force = false, xlab = L"A_1", ylab = L"A_2") 
end
