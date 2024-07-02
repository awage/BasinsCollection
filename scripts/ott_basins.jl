using DrWatson
@quickactivate 
using OrdinaryDiffEq:Vern9
using Attractors
using CairoMakie
# Equations of motion: E. Ott, et al. I Physica D 76 (1994) 384-410
function forced_particle!(du, u, p, t)
    γ=0.05  ; x̄ = 1.9  ; f₀=2.3  ; ω =3.5
    x₀=1. ; y₀=0.;
    x, y, dx, dy = u
    du[1] = dx
    du[2] = dy
    du[3] = -γ*dx -(-4*x*(1-x^2) + y^2) +  f₀*sin(ω*t)*x₀
    du[4] = -γ*dy -(2*y*(x+x̄)) +  f₀*sin(ω*t)*y₀
end


function _get_basins_ott(d)
    @unpack res = d
    xg = range(-2,2,length=res)
    yg = range(0.,2,length=res)
    diffeq = (reltol = 1e-9,  alg = Vern9())
    df = CoupledODEs(forced_particle!,rand(4),(0.0,20.0); diffeq)
    ω =3.5
    smap = StroboscopicMap(df, 2π/ω)
    psys = ProjectedDynamicalSystem(smap, [1,2], [0., 0,])
    mapper = AttractorsViaRecurrences(psys, (xg, yg); horizon_limit = 10)
    # The mapper search on a larger grid but we can focus on a tiny part of the 
    # phase space. (grid for recurrences and plotting are separate).   
    xg = range(0,1.2,length=res)
    yg = range(0.,1.2,length=res)
    bsn, att = basins_of_attraction(mapper, (xg, yg))
    grid = (xg,yg)
    return @strdict(bsn, xg, yg)
end


cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
# res = 1200 
params = @strdict res 
print_fig(params, "basins_riddle_ott", _get_basins_ott; cmap, xlab = L"x_0", ylab = L"y_0") 


