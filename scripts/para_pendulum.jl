using DrWatson
@quickactivate 
using Attractors
using OrdinaryDiffEq:Vern9
using CairoMakie
using LaTeXStrings
include(srcdir("print_fig.jl"))

# φ = (- b  ̇φ - c φ - m*L*(μ − A_2*ω^2 cos ωt) sin φ - m*L*A_1*ω^2 sin ωt cos φ)/(m*L^2) 
# Parametrically excited pendulum systems with several equilibrium positions. bifurcation analysis and rare attractors
# Alex v. Klokov∗ and Mikhail v. Zakrzhevsky
# http://dx.doi.org/10.1142/S0218127411030167
@inline @inbounds function pend_pforced(u, p, t)
    A_1, A_2, ω = p; x,y = u
    μ = 9.81; L = m = c = 1; b = 0.2
    du1 = y 
    du2 = -b*y -x -(μ - A_2*ω^2*cos(ω*t))*sin(x) - A_1*ω^2*sin(ω*t)*cos(x)
    return SVector{2}(du1, du2)
end


function compute_basins_pend_pforced(di::Dict)
    @unpack A_1, A_2, ω, res = di
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    ds = CoupledODEs(pend_pforced, rand(2), [A_1, A_2, ω]; diffeq)
    xg = yg = range(-14,14,length = 4001)
    smap = StroboscopicMap(ds, 2*pi/ω)
    xg = yg = range(-10,10,length = res)
    mapper = AttractorsViaRecurrences(smap, (xg, yg))
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


A_1 = 0.; A_2 = 4.1; ω = 1.5;
params = @strdict A_1 A_2 ω res
print_fig(params, "parametric_pendulum",compute_basins_pend_pforced; ylab= L"\dot{\theta}", xlab= L"\theta")
