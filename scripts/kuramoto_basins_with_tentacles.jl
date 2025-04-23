using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEqVerner
using ProgressMeter
using LinearAlgebra 
using Random

include(srcdir("print_fig.jl"))
include(srcdir("kuramoto_ring_def.jl"))


# Basins with Tentacles
# Yuanzhao Zhang and Steven H. Strogatz
# DOI: 10.1103/PhysRevLett.127.194101

# Generate a mapper for the Kuramoto oscillator, the mapper 
# takes the initial conditions on a subspace spanned by the 
# vector P1 and P2 such that  v = θ0 + α1 P1 + α2 P2, beeing 
# α1 and α2 two scalars. P1 and P2 and N dim vectors with 
# random components (shuffled 1s and 0s). This is done with 
# the function _get_ic. 
#  
# The mapper tracks the q-twisted states on a different
# projected space, it checks if the states converges to a 
# q-state with a lsq fit. 
function kuramoto_mapper(N) 
    ds = ring_kuramotos(N = N)
    P1 = shuffle([zeros(Int(N/2)) ; ones(Int(N/2))])
    P2 = shuffle([zeros(Int(N/2)) ; ones(Int(N/2))])
    θ₁ = 2π*10/N*(1:N)
    @show P1, P2, θ₁
    _get_ic(y) =   y[1].*P1 .+ y[2].*P2 .+ θ₁
    psys = ProjectedDynamicalSystem(ds, proj_fun, _get_ic)
    xg = range(-N/2,N/2; length = 1001)
    yg = range(-10,10; length = 11)
    mapper = AttractorsViaRecurrences(psys, (xg,xg); Δt = .1,   
    mx_chk_fnd_att = 1000,
    mx_chk_loc_att = 1000, Ttr = 150)
    return mapper
end

function compute_basin_tentacle(d::Dict) 
    @unpack N, res = d
    mapper = kuramoto_mapper(N)
    α1 = range(-π, π, length = res) 
    α2 = range(-π, π, length = res)
    bsn = @showprogress [mapper([a1, a2]) for a1 in α1, a2 in α2]
    grid = (α1,α2)
    att = extract_attractors(mapper)
    return @strdict(bsn, grid, res, att)
end

function kuramoto_featurizer(N) 
    ds = ring_kuramotos(N = N)
    m_q = round(N/4)
    templates = Dict( Int(u) => u for (k,u) in enumerate(range(-m_q, m_q, step = 1.)))
    mapper = AttractorsViaFeaturizing(ds, _proj_fun_f, GroupViaNearestFeature(templates); Ttr = 100., T = 200., Δt = 1.) 
    return mapper
end

function compute_basin_tentacle_feats(N, res)    
    mapper = kuramoto_featurizer(N)
    α1 = range(-π, π, length = res) 
    α2 = range(-π, π, length = res)
    P1 = shuffle([zeros(Int(N/2)) ; ones(Int(N/2))])
    P2 = shuffle([zeros(Int(N/2)) ; ones(Int(N/2))])
    θ₁ = 2π*10/N*(1:N)
    @show P1, P2, θ₁
    _get_ic(y) =   y[1].*P1 .+ y[2].*P2 .+ θ₁
    bsn = @showprogress [mapper(_get_ic([a1, a2])) for a1 in α1, a2 in α2]
    # ics = Vector{typeof(_get_ic([1.,1.]))}()
    # for x in α1, y in α2
    #     push!(ics, _get_ic([x,y]))
    # end
    # fs, labs = basins_fractions(mapper, StateSpaceSet(ics))
    # bas = reshape(labs,res,res)
    grid = (α1,α2)
    att = extract_attractors(mapper)
    # return @strdict(fs, labs, grid, res, att, bas)
    return @strdict(bsn, grid, res, att)
end

let res = 1200
    N = 40
    params = @strdict res N
    print_fig(params, "basins_tentacles", compute_basin_tentacle; force = false, xlab = L"\alpha_1", ylab = L"\alpha_2")
end

# dd = compute_basin_tentacle_feats(40, 100)

