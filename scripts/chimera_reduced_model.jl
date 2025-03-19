using DrWatson
@quickactivate 
using Attractors
using OrdinaryDiffEq:Vern9
using CairoMakie
using LaTeXStrings

# Basins of attraction for chimera states
# Erik A Martens1,2, Mark J Panaggio3,4 and Daniel M Abrams
# New J. Phys. 18 (2016) 022002
# doi:10.1088/1367-2630/18/2/022002
@inline @inbounds function chimera(u, p, t)
    ν, μ, β  = p
    s, d, ψ = u
    ρ1 = s + d
    ρ2 = s - d 
    dρ1 = (1-ρ1^2)/2*(μ*ρ1*sin(β) + ν*ρ2*sin(β-ψ))  
    dρ2 = (1-ρ2^2)/2*(μ*ρ2*sin(β) + ν*ρ1*sin(β+ψ))  
    ds = (dρ1 + dρ2)/2
    dd = (dρ1 - dρ2)/2
    dψ =  (1+ρ2^2)/(2*ρ2)*(μ*ρ2*cos(β) + ν*ρ1*cos(β+ψ)) - (1+ρ1^2)/(2*ρ1)*(μ*ρ1*cos(β) + ν*ρ2*cos(β-ψ))  
    return SVector{3}(ds, dd, dψ)
end

function affect!(integrator)
    uu = integrator.u
    if integrator.u[3] < 0
        set_state!(integrator, SVector(uu[1], uu[2], uu[3] + 2π))
        u_modified!(integrator, true)
    else
        set_state!(integrator, SVector(uu[1], uu[2], uu[3] - 2π))
        u_modified!(integrator, true)
    end
end

function compute_chimera(di::Dict)
    @unpack  res, ν, μ, β = di
    condition(u,t,integrator) = (integrator.u[3] < -π  || integrator.u[3] > π)
    cb = DiscreteCallback(condition,affect!)
    diffeq = (reltol = 1e-9,  alg = Vern9(), callback = cb)
    ds = CoupledODEs(chimera, rand(3), [ν, μ, β]; diffeq)
    xg = yg = range(-4,4, length = 10001)
    φ  = range(-5π,5π, length = 10001)
    mapper = AttractorsViaRecurrences(ds, (xg, yg, φ);
    consecutive_recurrences = 100,
    attractor_locate_steps = 1000)
    dg = range(-0.5, 0.5, length = res)
    ψg = range(0, 2π, length = res)
    grid = (dg,ψg)
    bsn = @showprogress [ mapper([y,0.56625, x]) for x in ψg, y in dg]
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid, res)
end


res = 1200; 
A = 0.1; β = 0.025
μ = (A+1)/2; ν = 1 - μ
params = @strdict res μ ν β 
print_fig(params, "chimera_states", compute_chimera; force = false, xlab = L"s", ylab = L"\psi") 
att = get_att(params, "chimera_states", compute_chimera) 
