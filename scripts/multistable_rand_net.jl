using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter
using LinearAlgebra 
using Random

include(srcdir("print_fig.jl"))

# PHYSICAL REVIEW E 90, 062710 (2014)
# Dynamics of random neural networks with bistable units
# M. Stern, 1,4,* H. Sompolinsky, 1,2,3 and L. F. Abbott 4,5
# DOI: 10.1103/PhysRevE.90.062710
function rand_net!(du, u, p, t)
    W,s,g,N = p
    f(x) = tanh(x)
    du .= -u .+ s*f.(u) .+ g/sqrt(N)*W*f.(u) 
    return nothing
end

function compute_basins_rand_net(di::Dict)
    @unpack  res, s, g, N = di
    rng = MersenneTwister(123);
    W = randn(rng,N,N)
    for k in 1:N; W[k,k] = 0.0; end;
   
    diffeq = (alg = Vern9(), reltol = 1e-6, maxiters = 1e8)
    ds = CoupledODEs(rand_net!, rand(N), (W,s,g,N); diffeq)
    yg =  collect(range(-20, 20; length = 5000))
    grid = ntuple(x -> yg, N)
    mapper = AttractorsViaRecurrences(ds, grid; 
    consecutive_recurrences = 100,
    Ttr = 4000)
    yr = range(-1., 1, length = res)
    xr = range(-1., 1, length = res)
    bsn = @showprogress [ mapper([x; y; zeros(N-2)]) for x in xr, y in yr]
    grid = (xr,yr)
    att = extract_attractors(mapper)
    return @strdict(bsn, grid, res, att)
end


# let 
res = 1200
# s = 4.5; g = 2.5; N = 100
s = 4.5; g = 0.4; N = 100
params = @strdict res s g N
print_fig(params, "basins_multi_rand_net", compute_basins_rand_net; force = false, xlab = L"x", ylab = L"y")
att = get_att(params, "basins_multi_rand_net", compute_basins_rand_net)
# end


    # s = 4.5; g = 0.4; N = 100
    # W = randn(N,N)
    # for k in 1:N; W[k,k] = 0.0; end;
    # diffeq = (alg = Vern9(), reltol = 1e-6, maxiters = 1e8)
    # ds = CoupledODEs(rand_net!, rand(N), (W,s,g,N); diffeq)
    # y,t = trajectory(ds, 100, [rand(2); zeros(N-2)])
