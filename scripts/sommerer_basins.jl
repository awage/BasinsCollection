using DrWatson
@quickactivate 
using OrdinaryDiffEq:Vern9
using Attractors
using CairoMakie
using ProgressMeter

# Sommerer, J. C. (1995). The end of classical determinism. Johns Hopkins APL Technical Digest, 16(4), 333.
function forced_particle!(du, u, p, t)
    γ = 0.632; f₀ = 1.0688  ; ω = 2.2136; 
    s = 20.; p = 0.098; k = 10.;
    x₀=1. ; y₀=0.;
    x, y, dx, dy = u
    du[1] = dx
    du[2] = dy
    du[3] = -γ*dx -(-4x*(1-x^2) + 2*s*x*y^2) +  f₀*sin(ω*t)*x₀
    du[4] = -γ*dy -(2*y*s*(x^2-p)+4*k*y^3) +  f₀*sin(ω*t)*y₀
end


function _get_basins_sommerer(d)
    @unpack res = d
    xg = range(-3,3,length = 3000)
    yg = range(-3,3,length = 3000)
    diffeq = (reltol = 1e-9,  alg = Vern9())
    df = CoupledODEs(forced_particle!,rand(4),(0.0,20.0); diffeq)
    ω = 2.2136
    smap = StroboscopicMap(df, 2π/ω)
    # psys = projected_integrator(smap, [1,2], [0., 0,])
    mapper = AttractorsViaRecurrences(smap, (xg, yg, xg, yg); 
            sparse = true)

    xg = range(-1.,1.,length=res)
    yg = range(-1.,1.,length=res)
    bsn = @showprogress [ mapper([x,y,0., 0.]) for x in xg, y in yg]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(bsn, grid, att)
end


res = 1000
params = @strdict res
print_fig(params, "basin_sommerer", _get_basins_sommerer)
