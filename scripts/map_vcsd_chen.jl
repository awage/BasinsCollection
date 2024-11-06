using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using ProgressMeter
include(srcdir("print_fig.jl"))

# Ostrovskii, V. Y., Rybin, V. G., Karimov, A. I., & Butusov, D. N. (2022). Inducing multistability in discrete chaotic systems using numerical integration with variable symmetry. Chaos, Solitons & Fractals, 165, 112794.
function map_chen!(dz, z, p, n)
xn, yn, zn = z
h = 0.01; a = 40; b = 3; c = 28 
S = p[1]
h1 = S*h
h2 = (1-S)*h
xns = (xn + h1*a*yn)/(1 + h1*a)
yns = (yn + h1*((c - a)*xns - xns*zn))/(1- h1*c)
zns = (zn + h1*xns*yns)/(1 + h1*b)
dzi = zns + h2*(xns*yns - b*zns)
dy = yns + h2*((c - a)*xns - xns*dzi + c*yns )
dx = xns + h2*a*(dy - xns) 
dz[1] = dx; dz[2] = dy; dz[3] = dzi
return
end

function compute_chen(di::Dict)
    @unpack S, res = di
    ds = DeterministicIteratedMap(map_chen!, rand(3), [S])
    yg = xg =  range(-100, 100., length = 1001)
    zg = range(-100, 100., length = 1001)
    # projection = [1,2]; complete_state = [20.]
    # pinteg = ProjectedDynamicalSystem(ds, projection, complete_state)
    # mapper = AttractorsViaRecurrences(pinteg, (xg, yg);
    mapper = AttractorsViaRecurrences(ds, (xg, yg, zg); 
    consecutive_attractor_steps = 10,
    consecutive_recurrences = 1000,
    attractor_locate_steps = 1000)
    yg = xg = range(-16., 16., length = res)
    grid = (xg, yg)
    bsn = @showprogress [ mapper([x,y,20.]) for x in xg, y in yg]
    # bsn = @showprogress [ mapper([x,y]) for x in xg, y in yg]
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid, res)
end

res = 1200; S = 0.77
params = @strdict res S
print_fig(params, "chen", compute_chen; xlab = L"x", ylab = L"y", force = false)
