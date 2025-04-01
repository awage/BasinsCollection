using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using ProgressMeter
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))

# title={Multistability in a dynamic Cournot game with three oligopolists},
# author={Agiza, Hamdy Nabih and Bischi, Gian Italo and Kopel, Michael},
# journal={Mathematics and Computers in Simulation},
# volume={51},
# number={1-2},
# pages={63--90},
# year={1999},
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

res = 1200; μ = 1.95; λ= 0.5
params = @strdict res μ λ 
cmap = ColorScheme([RGB(1,1,1), RGB(0.55,0.9,0.35), RGB(0.9,0.4,0.1),  RGB(0.50,0.24,1)] )
print_fig(params, "nash_equilibrium", compute_nash; xlab = L"q_1", ylab = L"q_2", force = false, cmap)
# att =  get_att(params, "nash_equilibrium", compute_nash; force = true)
