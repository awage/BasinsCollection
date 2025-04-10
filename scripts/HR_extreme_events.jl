using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using ProgressMeter
using OrdinaryDiffEq:Vern9
using ColorSchemes, Colors
include(srcdir("print_fig.jl"))

function hindmarshrose_two_coupled(u0=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
			a = 1.0, b = 3.0, c = 1.0, d = 5.0, r = 0.001, s = 4.0, xr = -1.6, I = 4.0,
			k1 = 0.05, k2 = 0.05, k_el = 0.0, xv = 2.0)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
	return CoupledODEs(hindmarshrose_coupled_rule, u0, [a, b, c, d, r, s, xr, I, k1, k2, k_el, xv]; diffeq)
end
function hindmarshrose_coupled_rule(u, p, t)
    function sigma(x)
        return @fastmath 1.0 / ( 1.0 + exp( -10.0 * ( x  - ( - 0.25 ) ) ) )
    end
    a, b, c, d, r, s, xr, I, k1, k2, k_el, xv = p
    x1, y1, z1, x2, y2, z2 = u

    du1 = y1 + b * x1 ^ 2 - a * x1 ^3 - z1 + I - k1 * ( x1 - xv ) * sigma(x2) + k_el * ( x2 - x1 )
    du2 = c - d * x1 ^2 - y1
    du3 = r * ( s * ( x1 - xr ) - z1 )

    du4 = y2 + b * x2 ^ 2 - a * x2 ^3 - z2 + I - k2 * ( x2 - xv ) * sigma(x1) + k_el * ( x1 - x2 )
    du5 = c - d * x2 ^2 - y2
    du6 = r * ( s * ( x2 - xr ) - z2 )
    return SVector(du1, du2, du3, du4, du5, du6)
end

function compute_HR(di::Dict)
    @unpack res = di
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    ds = hindmarshrose_two_coupled()
    yg = range(-30, 30; length = 1001)
    grid = ntuple(x -> yg, 6)
    mapper = AttractorsViaRecurrences(ds, grid; Δt = 1., maximum_iterations = Int(1e9),
    consecutive_recurrences = 12000,
    Ttr = 1000,
    attractor_locate_steps = 1000,
    consecutive_attractor_steps = 25)
    y1r = range(-4, 4, length = res)
    y2r = range(-4, 4, length = res)
    bsn = @showprogress [ mapper([y1, 1, 5, y2, 1, 5]) for y1 in y1r, y2 in y2r] 
    grid = (y1r,y2r); att = mapper.bsn_nfo.attractors
    return @strdict(bsn, grid, att, res)
end

res = 400
params = @strdict res 
cmap = ColorScheme([ RGB(0.9,0.2,0.1), RGB(1,1,1) ] )
print_fig(params, "coupled_HR", compute_HR; force = true, xlab = L"z_1", ylab = L"z_2", cmap) 
att = get_att(params, "coupled_HR", compute_HR; force = false) 


    ds = hindmarshrose_two_coupled()
