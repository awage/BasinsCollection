using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors
using StaticArrays
using ProgressMeter
using OrdinaryDiffEqVerner
include(srcdir("print_fig.jl"))

# Multistability, phase diagrams, and intransitivity in the Lorenz-84 low-order atmospheric circulation model
# Chaos 18, 033121 (2008); https://doi.org/10.1063/1.2953589
@inline @inbounds function lorenz84(u, p, t)
    F = p[1]; G = p[2]; a = p[3]; b = p[4];
	x = u[1]; y = u[2]; z = u[3];
    dx = -y^2 -z^2 -a*x + a*F
    dy = x*y - y - b*x*z +G
	dz = b*x*y + x*z -z
    return SVector{3}(dx, dy, dz)
end

function compute_lorenz84(di)
    @unpack F, G, a, b, res = di
    diffeq = (alg = Vern9(), reltol = 1e-7, maxiters = 1e8)
    ds =  CoupledODEs(lorenz84, rand(3), [F,G,a,b]; diffeq)
    xg = yg = zg = range(-4, 4,length = 250)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); 
                consecutive_recurrences = 2000, 
                attractor_locate_steps = 2000,
                consecutive_attractor_steps = 100)
    xg = range(0.5,1.5,length = res)
    yg = range(-1.5,-0.5,length = res)
    bsn = @showprogress [ mapper([x,y,0.]) for x in xg, y in yg]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(bsn, att, grid, F, G, a, b, res)
end

let
    F=6.846; G=1.287; a=0.25; b=4.; res = 1200
    params = @strdict F G a b res
    print_fig(params, "lorenz84", compute_lorenz84; ylab = L"y", xlab = L"x", force = false)
end
