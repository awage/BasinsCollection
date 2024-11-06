using DrWatson
@quickactivate 
using Attractors
using OrdinaryDiffEq
using LaTeXStrings
using CairoMakie
include(srcdir("print_fig.jl"))

# Adaptive synapse-based neuron model with
# heterogeneous multistability and riddled basins
#  Chaos 32, 123101 (2022); doi: 10.1063/5.0125611
# H. Bao, J. Zhang,1 N. Wang, N. V. Kuznetsov
# and B. C. Bao1,a)
function neuron_synapse!(du, u, p, t)
    x,y = u; c, B, g = p
    H(x) = B*sin(g*x) 
    du[1] = -x + H(x)*H(y)  + sin(2π*t)
    du[2] = -c*y + c*H(x)^2
end


function compute_neuron_synapse(di::Dict)
    @unpack  c, B, g, res = di
    diffeq = (alg = Vern9(), reltol = 1e-6, maxiters = 1e8)
    ds = CoupledODEs(neuron_synapse!, zeros(2), [c, B, g] ; diffeq)
    pstrob = StroboscopicMap(ds, 2π)
    y1 = range(-15, 15, length = 1000)
    y2 = range(-15, 15, length = 1000)
    mapper = AttractorsViaRecurrences(pstrob, (y1,y2))
        # consecutive_recurrences = 200) 
    y1 = range(-10, 10, length = res)
    y2 = range(-10, 10, length = res)
    bsn, att = basins_of_attraction(mapper, (y1,y2))
    grid = (y1,y2)
    return @strdict(bsn, att, grid,  res)
end

c = 1.8; B = 2; g = 1.7
res = 1200
params = @strdict c B g res
print_fig(params, "neuron_synapse", compute_neuron_synapse; ylab = L"y", xlab = L"x", force = true)


