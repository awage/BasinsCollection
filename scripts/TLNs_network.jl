using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter
using LinearAlgebra 

include(srcdir("print_fig.jl"))
# Core motifs predict dynamic attractors in combinatorial threshold-linear networks
#     Caitlyn Parmelee,  Samantha Moore,Katherine Morrison , Carina Curto
#  https://doi.org/10.1371/journal.pone.0264456 
# Basins are not computed in the paper but there is 
# multistability. 
function TLN_net!(du, u, p, t)
    W,b = p
    du .= -u .+ max.(W*u .+ b, 0)
    return nothing
end




function compute_basins_tlns(di::Dict)
    @unpack  res,θ = di
    ε = 0.25; δ = 0.5
    W = ones(9,9) - diagm(ones(9))
    W = W .* (-1 - δ)
    W[2,1] = W[4,1] = W[8,1] = W[9,1] = -1 + ε
    W[6,2] = W[5,2] = -1 + ε
    W[2,3] = W[4,3] = -1 + ε
    W[5,4] = W[8,4] = -1 + ε
    W[1,5] = W[3,5] = W[6,5] = -1 + ε
    W[3,6] =  -1 + ε
    W[8,7] = W[1,7] = -1 + ε
    W[7,8] = W[1,8] = W[9,8] = W[4,8] = -1 + ε
    W[1,9] = W[2,9] = W[8,9] =  -1 + ε
    b = ones(9,1).*θ
    diffeq = (alg = Vern9(), reltol = 1e-6, maxiters = 1e8)
    ds = CoupledODEs(TLN_net!, rand(9), [W,b]; diffeq)
    yg =  collect(range(0, 10; length = 5000))
    grid = ntuple(x -> yg, 9)
    mapper = AttractorsViaRecurrences(ds, grid; 
    attractor_locate_steps = 1000)
    yr = range(0., 4, length = res)
    xr = range(0., 4, length = res)
    bsn = @showprogress [ mapper([zeros(5); x;y; 0.; 0.]) for x in xr, y in yr]
    # bsn = @showprogress [ mapper([x;y;zeros(7)]) for x in xr, y in yr]
    grid = (xr,yr)
    att = extract_attractors(mapper)
    return @strdict(bsn, grid, res, att)
end


let res = 1200
θ = 1.; 
params = @strdict res θ
print_fig(params, "basins_tlns", compute_basins_tlns; force = false)
end

