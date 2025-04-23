using DrWatson
@quickactivate # exports Attractors, GLMakie and other goodies in `src`
using Attractors
using OrdinaryDiffEqVerner
using CairoMakie
using LaTeXStrings
using Colors,ColorSchemes
using ProgressMeter
include(srcdir("print_fig.jl"))

# Miriam Leon, Mae L. Woods, Alex J. H. Fedorec, Chris P. Barnes
# "A computational method for the investigation of multistable systems and its application to genetic switches"
# [10.1186/s12918-016-0375-z](https://doi.org/10.1186/s12918-016-0375-z)

function ludp(u, p, t)
    # d = p[1]; F = p[2]; omega = p[3]
    xx = 250
    λxx = 10
    nxx = 4
    Hxx(x) = (1-λxx)/(1 + (x/xx)^nxx) + λxx
    Hyy(x) = Hxx(x)

    xy = 500
    λxy = 0.5
    nxy = 2
    Hxy(x) =  (1-λxy)/(1 + (x/xy)^nxy) + λxy

    yx = 500
    λyx = 0.1
    nyx = 2
    Hyx(x) =  (1-λyx)/(1 + (x/yx)^nyx) + λyx

    gx = 40; gy = 40; kx = 0.55; ky = 0.55

    x,y = u 
    du1 = gx*Hxy(y)*Hxx(x) - kx*x
    du2 = gy*Hyx(x)*Hyy(y) - ky*y
    return SVector{2}(du1, du2)
end



function compute_basins_ludp(di::Dict)
    @unpack res = di
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    ds = CoupledODEs(ludp, rand(2) ; diffeq)
    xg = yg = range(0,2000,length = 400)
    mapper = AttractorsViaRecurrences(ds, (xg, yg))
    xg = yg = range(0,350,length = res)
    grid = (xg, yg)
    bsn = @showprogress [ mapper([x,y]) for x in xg, y in yg]
    att = mapper.bsn_nfo.attractors
    return @strdict(bsn, att, grid,  res)
end



res = 400
    params = @strdict res
    print_fig(params, "ludp", compute_basins_ludp; ylab = L"y", xlab = L"x", force = false) 
