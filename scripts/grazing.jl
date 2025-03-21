using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes
# Physica D 130 (1999) 43–57
# Grazing bifurcations and basins of attraction in an impact-friction
# oscillator
# L.N. Virgin ∗, C.J. Begley
function g(x) 
    σl = -1.5; σr = 1.5; ρ = 1.;
    if x ≥ σr 
        return(ρ^2*(x - σr))
    elseif σl ≤ abs(x) ≤ σr
        return(0)
    elseif x ≤ σl
        return(ρ^2*(x + σl))
    end
end

# Equations of motion:
function grazing(u, p, t)
    h = 0.05 ; γ = 1.8; α = 1.; ff = 0.01
    f(x) = α*ff*sign(x)
    du1 = u[2]
    du2 = - 2*h*u[2] - f(u[2]) -g(u[1]) + 2*h*γ*cos(γ*t) + sin(γ*t)
    return SVector{2}(du1, du2)
end



function compute_grazing(di::Dict)
    @unpack res = di
    diffeq = (reltol = 1e-8,  alg = Vern9(),)
    df = CoupledODEs(grazing, rand(2); diffeq)
    xg = range(-4.,4.,length = 50001)
    yg = range(-4.,4.,length = 50001)
    mapper = AttractorsViaRecurrences(df, (xg, yg);
            consecutive_attractor_steps = 4,
            # consecutive_basin_steps = 400,
            consecutive_recurrences = 4000,
            attractor_locate_steps = 2000)
    xg = range(-1,1,length = res)
    yg = range(-2.,2.,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


res = 10
params = @strdict res
print_fig(params, "grazing", compute_grazing; ylab= L"\dot{\theta}", xlab= L"\theta", force = true)
att = get_att(params, "grazing", compute_grazing)

