using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))


# @article{lorenz1989computational,
#   title={Computational chaos-a prelude to computational instability},
#   author={Lorenz, Edward N},
#   journal={Physica D: Nonlinear Phenomena},
#   volume={35},
#   number={3},
#   pages={299--317},
#   year={1989},
#   publisher={Elsevier}
# }
# Equations of motion:
function lorenz_comp(u, p, t)
    τ = 1.5; a = 0.36
    x,y = u
    du1 = (1 + a*τ)*x - τ*x*y
    du2 = (1-τ)*y + τ*x^2
    return SVector{2}(du1, du2)
end



function compute_lorenz_comp(di::Dict)
    @unpack res = di
    df = DeterministicIteratedMap(lorenz_comp, rand(2))
    xg = range(-4.,4.,length = 501)
    yg = range(-4.,4.,length = 501)
    mapper = AttractorsViaRecurrences(df, (xg, yg);
            consecutive_attractor_steps = 4,
            # consecutive_basin_steps = 400,
            consecutive_recurrences = 500,
            attractor_locate_steps = 200)
    xg = range(-2,2,length = res)
    yg = range(-3.,3.,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


res = 1200
params = @strdict res
cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(1,0.46, 0.46), RGB(0.34,0.34,1), RGB(0.1,0.1,0.1) ] )
print_fig(params, "lorenz_complexity", compute_lorenz_comp;  force = false, cmap)
att = get_att(params, "lorenz_complexity", compute_lorenz_comp)

