using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEqVerner
using ProgressMeter
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))
# This system is problematic: very very long transient (t > 2000 sometimes) that can be mistaken with attractors.

# ```
# Rikitake's dynamo [^Rikitake1958] is a system that tries to model the magnetic
# reversal events by means of a double-disk dynamo system.

# [^Rikitake1958]: T. Rikitake Math. Proc. Camb. Phil. Soc. **54**, pp 89–105, (1958)
# """
function rikitake_rule(u, p, t)
    μ, α = p
    x,y,z = u
    xdot = -μ*x + y*z
    ydot = -μ*y + x*(z - α)
    zdot = 1 - x*y
    return SVector{3}(xdot, ydot, zdot)
end


function compute_rikitake(di::Dict)
    @unpack μ, α, res = di
    diffeq = (alg = Vern9(), reltol = 1e-6, maxiters = 1e8)
    ds = CoupledODEs(rikitake_rule, rand(3), [μ, α]; diffeq)
    xg = yg = range(-5,5,length=1001)
    psys = ProjectedDynamicalSystem(ds, [1,2], [0.])
    mapper = AttractorsViaRecurrences(psys, (xg, yg) ; Δt = 0.1,
        mx_chk_fnd_att = 50000,
        mx_chk_loc_att = 10000)
    xg = yg = range(-2.5, 2.5, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg))
    att = mapper.bsn_nfo.attractors
    grid = (xg,yg)
    return @strdict(bsn, att, grid, μ, α, res)
end

μ = 0.5; α = 1.;
params = @strdict res μ α
cmap = ColorScheme([RGB(1,1,1), RGB(0.9,0.15,0.15),RGB(0.95,0.95,0.95), RGB(0.95,0.95,0.95),   RGB(0.95,0.95,0.95),RGB(0.95,0.95,0.95),  RGB(0.15,0.9,0.15), RGB(0.9,0.4,0.1),  RGB(1.0,1.0,1.)] )
print_fig(params, "rikitake", compute_rikitake; ylab = L"y", xlab = L"x", force = false, cmap)
