using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using ProgressMeter
include(srcdir("print_fig.jl"))

# @article{do2008multistability,
#   title={Multistability and arithmetically period-adding bifurcations in piecewise smooth dynamical systems},
#   author={Do, Younghae and Lai, Ying-Cheng},
#   journal={Chaos: An Interdisciplinary Journal of Nonlinear Science},
#   volume={18},
#   number={4},
#   year={2008},
#   publisher={AIP Publishing}
# }
function map_lai!(dz, z, p, n)
    x,y = z
    a = -2; b = -0.95;
    c = -b/a; d = b; μ = -1.0
    if x ≤ 0 
        dx = a*x + y + μ
        dy = b*x
    else 
        dx = c*x + y + μ
        dy = d*x
    end
    dz[1] = dx; dz[2] = dy; 
    return
end

function compute_lai(di::Dict)
    @unpack  res = di
    ds = DeterministicIteratedMap(map_lai!, rand(2))
    yg = xg =  range(-2, 2., length = 5001)
    mapper = AttractorsViaRecurrences(ds, (xg, yg); 
            consecutive_attractor_steps = 2,
            consecutive_basin_steps = 400,
            consecutive_recurrences = 1000)
            # attractor_locate_steps = 2000)
        # consecutive_recurrences = 1000)
    yg = xg = range(-1.5, 2.5, length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, res)
end

res = 1200; 
params = @strdict res 
cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(0.1,0.1,0.1), RGB(1,0.46, 0.46), RGB(0.34,0.34,1)] )
print_fig(params, "lai_pcw", compute_lai; xlab = L"x", ylab = L"y", cmap, force = false)
# att =  get_att(params, "lai_pcw", compute_lai)
