using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes
include(srcdir("print_fig.jl"))

# @article{di1997secondary,
#   title={Secondary bifurcations and high periodic orbits in voltage controlled buck converter},
#   author={Di Bernardo, Mario and Fossas, Enric and Olivar, Gerard and Vasca, Francesco},
#   journal={International Journal of Bifurcation and Chaos},
#   volume={7},
#   number={12},
#   pages={2755--2771},
#   year={1997},
#   publisher={World Scientific}
# }
# Equations of motion:
function buck(u, p, t)
    R = 22; C = 47e-6; L = 20e-3; γ = 11.75
    η = 1309.52; T= 400e-6; Vin = p[1]
    vr = γ + η*mod(t,T)  
    ur = u[1] > vr ? 0. : 1.  
    du1 = -1/(R*C)*u[1] + 1/C*u[2]
    du2 = -1/L*u[1] + Vin/L*ur
    return SVector{2}(du1, du2)
end


function compute_buck(di::Dict)
    @unpack res, Vin = di
    T = 400e-6
    diffeq = (reltol = 0, abstol= 1e-14, alg = Vern9())
    df = CoupledODEs(buck, rand(2), [Vin]; diffeq)
    smap = StroboscopicMap(df, T) 
    xg = range(-20,20.0,length = 5001)
    yg = range(-4.,4.,length = 5001)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    xg = range(10,14,length = res)
    yg = range(0.5,1,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end

res = 600; Vin = 30.10
params = @strdict res Vin
print_fig(params, "buck", compute_buck; ylab= L"V", xlab= L"i", force = true)


