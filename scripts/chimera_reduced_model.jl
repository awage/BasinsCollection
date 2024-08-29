using DrWatson
@quickactivate 
using Attractors
using OrdinaryDiffEq:Vern9
using CairoMakie
using LaTeXStrings

@inline @inbounds function chimera(u, p, t)
    ν, μ, β  = p
    s, d, ψ = u
    # ρ1, ρ2, ψ = u
    ρ1 = s + d
    ρ2 = s - d 
    dρ1 = (1-ρ1^2)/2*(μ*ρ1*sin(β) + ν*ρ2*sin(β-ψ))  
    dρ2 = (1-ρ2^2)/2*(μ*ρ2*sin(β) + ν*ρ1*sin(β+ψ))  
    dψ =  (1+ρ2^2)/(2*ρ2)*(μ*ρ2*cos(β) + ν*ρ1*cos(β+ψ)) - (1+ρ1^2)/(2*ρ1)*(μ*ρ1*cos(β) + ν*ρ2*cos(β-ψ))  
    ds = (dρ1 + dρ2)/2
    dd = (dρ1 - dρ2)/2
    return SVector{3}(ds, dd, dψ)
    # return SVector{3}(dρ1, dρ2, dψ)
end


function compute_chimera(di::Dict)
    @unpack  res, ν, μ, β = di
    diffeq = (;reltol = 1e-7, alg = Vern9(), maxiters = 1e6)
    ds = CoupledODEs(chimera, rand(3), [ν, μ, β]; diffeq)
    xg = yg = range(-2,2, length = 10000)
    φ  = range(-2π,2π, length = 10000)
    mapper = AttractorsViaRecurrences(ds, (xg, yg, φ);
    consecutive_recurrences = 1000,
    attractor_locate_steps = 1000)
    dg = range(-0.5, 0.5, length = res)
    ψg = range(-π, π, length = res)
    grid = (dg,ψg)
    # bsn = @showprogress [ mapper([x,1-(μ-ν), y]) for x in dg, y in ψg]
    bsn = @showprogress [ mapper([y,0.56625, x]) for x in ψg, y in dg]
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid, res)
end


res = 1200; 
A = 0.1; β = 0.025
μ = (A+1)/2; ν = 1 - μ
params = @strdict res μ ν β 
print_fig(params, "chimera_states", compute_chimera; force = true) 
