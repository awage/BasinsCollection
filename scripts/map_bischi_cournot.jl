using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
include(srcdir("print_fig.jl"))

# @article{bischi2000multistability,
#   title={Multistability and cyclic attractors in duopoly games},
#   author={Bischi, Gian Italo and Mammana, Cristiana and Gardini, Laura},
#   journal={Chaos, Solitons \& Fractals},
#   volume={11},
#   number={4},
#   pages={543--564},
#   year={2000},
#   publisher={Elsevier}
# }
function cournot_2d!(dz, z, p, n)
    xn = z[1]; yn = z[2]
    a, ε = p
    dz[1] = μ1*yn*(1-yn)
    dz[2] = μ2*xn*(1-xn)
    return
end


function compute_cournot_2d(di::Dict)
    @unpack μ1, μ2, res = di
    ds = DeterministicIteratedMap(cournot_2d!, [1.0, 0.0], [μ1, μ2])
    yg = xg = range(-2., 2., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg);     
        consecutive_recurrences = 1000)
    yg = xg = range(0.0, 1., length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end


res = 1200
μ1 = 3.55; μ2 = 3.55
params = @strdict res μ1 μ2
# cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(0.34,0.34,1), RGB(1,0.46,0.46), RGB(0.1,0.1,0.1) ] )
print_fig(params, "cournot2d", compute_cournot_2d; force = true) 
# att = get_att(params, "cournot2d", compute_cplog) 

