using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors
using ProgressMeter
using OrdinaryDiffEq:Vern9
include(srcdir("print_fig.jl"))
function biorhythm!(du, u, p, t)
    σ1 = 10; σ2 = 10; L_1 = 5e8; L_2 = 100
    q1 =50; q2 = 0.02; d = 1e-6
    ks = p[1]; v_Km1 = 0.45 
    α, β, γ = u
    Φ = α*(1 + α)*(1 + β)^2/(L_1 + (1 + α)^2*(1 + β)^2)
    η = β*(1 + d*β)*(1 + γ)^2/(L_2 + (1 + d*β)^2*(1 + γ)^2)
    du[1] = v_Km1 - σ1*Φ
    du[2] = q1*σ1*Φ - σ2*η
    du[3] = q2*σ2*η - ks*γ
    return nothing
end

function compute_biorhythm(di)
    @unpack ks, res = di
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    ds =  CoupledODEs(biorhythm!, rand(3), [ks]; diffeq)
    # xg = yg = zg = range(-1, 1000,length = 5000)
    pow = 2; xg = range(0, 2000^(1/pow); length = 5001).^pow
    mapper = AttractorsViaRecurrences(ds, (xg,xg,xg);
                consecutive_recurrences = 2001,
                attractor_locate_steps = 2000)
                # consecutive_attractor_steps = 100)
    xg = range(0.,250,length = res)
    yg = range(0.,100,length = res)
    bsn = @showprogress [ mapper([x,y,1.]) for x in xg, y in yg]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(bsn, att, grid,  res)
end

res = 1200; 

# for ks in 2.0:0.05:2.5
ks = 1.99
params = @strdict  ks res
print_fig(params, "biorhythm", compute_biorhythm; ylab = L"y", xlab = L"x", force = false)
att = get_att(params, "biorhythm", compute_biorhythm)
# end


# ks =2 
# diffeq = (alg = Vern9(), reltol = 1e-7, maxiters = 1e8)
# ds =  CoupledODEs(biorhythm!, 100*rand(3), [ks]; diffeq)
# y,t = trajectory(ds, 10000, [rand(2)*101; 1])
# lines(t,y[:,1])
