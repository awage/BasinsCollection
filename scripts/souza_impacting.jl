using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes

function H(x) 
    gg = 0.63
    if x ≥ gg 
        return 1
    elseif x < gg
        return 0
    end
end

# Equations of motion:
function impact_oscillator(u, p, t)
    @inbounds begin
    ξ = 0.02 ; β = 0.5; ω = 0.417893; α = 20; gg = 0.63 
    f(y) = y
    du1 = u[2]
    du2 = -2*ξ*u[2] -f(u[1]) -α*H(u[1])*(u[1]-gg) + β*cos(ω*t)
    return SVector{2}(du1, du2)
    end
end


function compute_impact(di::Dict)
    @unpack res = di
    ω = 0.417893;
    diffeq = (reltol = 1e-8,  alg = Vern9(),)
    df = CoupledODEs(impact_oscillator, rand(2); diffeq)
    xg = range(-2.,2.,length = 5001)
    yg = range(-2.,2.,length = 5001)
    smap = StroboscopicMap(df, 2*pi/ω)
    mapper = AttractorsViaRecurrences(smap, (xg, yg);
            consecutive_attractor_steps = 3,
            consecutive_basin_steps = 400,
            consecutive_recurrences = 1000)
            # attractor_locate_steps = 500)
    xg = range(-1,1,length = res)
    yg = range(-1.,1.,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


res = 1200
params = @strdict res
print_fig(params, "impacting", compute_impact; ylab= L"\dot{y}", xlab= L"y", force = false)
# att = get_att(params, "gear_rattle", compute_gear_rattle)

