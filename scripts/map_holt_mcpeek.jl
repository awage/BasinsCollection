using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using ProgressMeter
include(srcdir("print_fig.jl"))
     # Chaos, Solitons and Fractals 12 (2001) 301±311
# Dynamics with riddled basins of attraction in models of
# interacting populations
# Bernard Cazelles 
function map_cazelle!(dz, z, p, n)
x1, y1, x2, y2 = z
dx = 0.175; dy = 0.002; k1 = 100.0; k2 = 75.0; r1 = 3.0; r2 = 3.5
dz[1] = (1-dx)*x1*exp(r1*(1-(x1+y1)/k1)) + dx*x2*exp(r2*(1-(x2+y2)/k2))
dz[2] = (1-dy)*y1*exp(r1*(1-(x1+y1)/k1)) + dy*y2*exp(r2*(1-(x2+y2)/k2))
dz[3] = (1-dx)*x2*exp(r2*(1-(x2+y2)/k2)) + dx*x1*exp(r1*(1-(x1+y1)/k1))
dz[4] = (1-dy)*y2*exp(r2*(1-(x2+y2)/k2)) + dy*y1*exp(r1*(1-(x1+y1)/k1))
# dz = dx1, dy1, dx2, dy2
return nothing
end

function compute_cazelle(di::Dict)
    @unpack res = di
    ds = DeterministicIteratedMap(map_cazelle!, rand(4))
    y1 = x1 =  range(-1, 10., length = 1001)
    y2 = x2 =  range(-1, 10., length = 1001)
    mapper = AttractorsViaRecurrences(ds, (x1,x2,y1,y2); 
    consecutive_attractor_steps = 10,
    consecutive_recurrences = 1000,
    attractor_locate_steps = 1000)
    yg = xg = range(0, 2., length = res)
    grid = (xg, yg)
    bsn = @showprogress [ mapper([x,y,x/2,y/2]) for x in xg, y in yg]
    # bsn = @showprogress [ mapper([x,y]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid, res)
end

res = 1200; 
params = @strdict res 
cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
print_fig(params, "cazelle", compute_cazelle; xlab = L"x", ylab = L"y", force = false, cmap)
