using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors
using StaticArrays
using ProgressMeter
using OrdinaryDiffEqVerner
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))

# https://doi.org/10.3390/math10111914
@inline @inbounds function gyrostat(u, p, t)
    b11 = 2; b12 = 0.7933; b13 = 0.1914
    b21 = 1.19; b22 = 3.48; 
    b31 = 0.5742; b33 = 5.8
    F1m = 1/3; F2m = -1; F3m = 1; 
    L1m =0; L2m = 0; L3m = 22.8 
	x1 = u[1]; x2 = u[2]; x3 = u[3];
    dx1 = -b11*x1 -b12*x2 +b13*x3 + F1m*x2*x3 + L1m
    dx2 = b21*x1 + b22*x2 + F2m*x1*x3 + L2m
	dx3 = -b31*x1 - b33*x3 + F3m*x1*x2 + L3m
    return SVector{3}(dx1, dx2, dx3)
end

function compute_gyrostat(di)
    @unpack  res = di
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    ds =  CoupledODEs(gyrostat, rand(3); diffeq)
    xg = yg = zg = range(-30, 30,length = 1501)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); Ttr = 2000,
                                      consecutive_recurrences = 1000)
    xg = range(-50,50,length = res)
    yg = range(-50,50,length = res)
    bsn = @showprogress [ mapper([x,y,0.]) for x in xg, y in yg]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(bsn, att, grid,  res)
end

res = 1200
params = @strdict res
cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.25,0.2)] )
print_fig(params, "gyrostat", compute_gyrostat; ylab = L"x_2", xlab = L"x_1", force = false, cmap)
att = get_att(params, "gyrostat", compute_gyrostat)

