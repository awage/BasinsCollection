using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))

# Bashkirtseva, I.A.; Ryashko, L.B.; Pisarchik, A.N. 
# Dragon Intermittency at the Transition to Synchronization in Coupled Rulkov
# Neurons. Mathematics 2025, 13, 415.
# https://doi.org/10.3390/math13030415
function cplog_rulkov(dz, z, p, n)
    xn = z[1]; yn = z[2]
    α, γ1, γ2, σ = p
    f(γ, x) = α/(1+x^2) + γ 
    dz[1] = f(γ1, xn) + σ*(yn - xn) 
    dz[2] =  f(γ2, yn) + σ*(xn-yn)   
    return
end


function compute_cpld_rulkov(di::Dict)
    @unpack γ1, γ2, α, σ, res = di
    ds = DeterministicIteratedMap(cplog_rulkov, [1.0, 0.0], [α, γ1, γ2, σ])
    yg = xg = range(-3., 3., length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,yg);     
        consecutive_recurrences = 1000, Ttr = 100)
    yg = xg = range(-2, 2, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end


res = 1200
Δ = 0.0; γ1 = -1.75; γ2 = γ1 + Δ; α = 4.1 ; σ = 0.035
params = @strdict res γ1 γ2  α Δ σ
cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(0.34,0.34,1), RGB(1,0.46,0.46), RGB(0.1,0.1,0.1) ] )
print_fig(params, "coupled_rulkov", compute_cpld_rulkov; force = false, cmap) 


