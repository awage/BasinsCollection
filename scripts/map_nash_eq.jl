using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using ProgressMeter
include(srcdir("print_fig.jl"))

# International Journal of Bifurcation and Chaos, Vol. 19, No. 1 (2009) 203–224
# c World Scientific Publishing Company
# BASINS OF ATTRACTION IN NONSMOOTH
# MODELS OF GEAR RATTLE
# JOANNA F. MASON∗
# MACSI, Department of Mathematics & Statistics,
# PETRI T. PIIROINEN
# R. EDDIE WILSON and MARTIN E. HOMER
function nash_map!(dz, z, p, n)
    q1, q2, q3 = z
    μ,λ = p
    dz[1] = (1 - λ)*q1 + λ*μ*(q2*(1-q2)+q3*(1-q3))
    dz[2] = (1 - λ)*q2 + λ*μ*(q3*(1-q3)+q1*(1-q1))
    dz[3] = (1 - λ)*q3 + λ*μ*(q1*(1-q1)+q2*(1-q2))
    return
end

function compute_nash(di::Dict)
    @unpack μ, λ, res = di
    ds = DeterministicIteratedMap(nash_map!, rand(3), [μ,λ])
    yg = xg = zg =  range(0.5, 1., length = 45001)
    mapper = AttractorsViaRecurrences(ds, (xg, yg, zg))
    yg = xg = range(-0.5, 1.8, length = res)
    bsn = @showprogress [ mapper([x,y,1-1/(2*μ)]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end

res = 200; μ = 1.95; λ= 0.5
params = @strdict res μ λ 
print_fig(params, "nash_equilibrium", compute_nash; xlab = L"q_1", ylab = L"q_2", force = true)
# att =  get_att(params, "nash_equilibrium", compute_nash; force = true)
