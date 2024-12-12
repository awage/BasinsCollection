using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
# using OrdinaryDiffEq:Vern9
using OrdinaryDiffEq
using ProgressMeter
include(srcdir("print_fig.jl"))

# Generalized multistability and its control in a laser featured
# Riccardo Meucci, Jean Marc Ginoux, Mahtab Mehrabbeik, Sajad Jafari, Julien Clinton Sprott
# Chaos 32, 083111 (2022)
# https://doi.org/10.1063/5.0093727
function CO2_laser!(du, u, p, t)
    m,f = p
    k = 12; B0 = 0.05; γ = 0.0025; 
    α = 0.002; p0 = 1.252 
    x,y = u
    # du[1] = -x*(1+k*(B0+m*sin(2π*f*t))^2 - y)
    # du[2] = -γ*y -α*x*y + γ*p0

    # Change of variables v = log(x) to avoid 
    # numerical instabilities   
    du[1] = -(1+k*(B0+m*sin(2π*f*t))^2 - y)
    du[2] = -γ*y -α*exp(x)*y + γ*p0
    return nothing
end


function compute_basins_CO2_lasers(di::Dict)
    @unpack  res, m ,f = di
    diffeq = (; reltol = 1e-6, abstol = 1e-8, alg = Vern9(), maxiters = 1e6)
    ds = CoupledODEs(CO2_laser!, rand(2), [m,f]; diffeq)
    smap = StroboscopicMap(ds, 1/f)
    yg = xg =  range(-30, 30; length = 20000)
    # pow = 3; xg = range(-20, 20^(1/pow); length = 5000).^pow
    # yg =  collect(range(-20, 20; length = 5000))
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(smap, grid; 
    consecutive_recurrences = 2000)
    # Δt = 1/(10f), force_non_adaptive = true,
    # consecutive_recurrences = 2000)
    yr = range(0.1, 4, length = res)
    # xr = log.(collect(range(0.1, 15, length = res)))
    xr = range(0.1, 8, length = res)
    bsn = @showprogress [ mapper([x,y]) for x in xr, y in yr]
    # xr = range(0.1, 15, length = res)
    grid = (xr,yr)
    att = extract_attractors(mapper)
    return @strdict(bsn, grid, res, att)
end


let res = 1200
    m = 0.02; f = 0.005
    params = @strdict res m f
    cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
    print_fig(params, "basins_CO2_laser", compute_basins_CO2_lasers; force = false, xlab = L"\log x", ylab = L"y", cmap)
end

