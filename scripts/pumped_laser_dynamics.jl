using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter
include(srcdir("print_fig.jl"))
# using Plots


# Control of basins of attraction in a multistable fiber laser
# A.N. Pisarchik a,∗ , R. Jaimes-Reategui b
# Physics Letters A 374 (2009) 228–234
# doi:10.1016/j.physleta.2009.10.061

function pumped_laser!(du, u, p, t)
    L = 70; Tr = 8.7; r0 = 1.5e-4; σ12 = 3e-21; 
    τ = 1e-2; ξ1 = 2; ξ2= 0.4; N0 = 1.3e20; R = 0.8
    w0 = 3.5e-4; γ0 = 0.038; λg = 1.56e-4; β = 0.5
    # α0 = N0*σ12 
    α0 = 0.4
    # αth = γ0 + 1/(2*L)*log(1/R)
    αth = 3.92e-2
    # rw = 1 - exp(-2*(r0/w0)^2)
    rw = 0.308
    P,N = u
    md = 0.8; fd = 70.2e3
    Pp = 7.4e19*(1-md*sin(2*π*fd*t))
    Psp = N*1e-3/(τ*Tr)*(λg/w0)^2*r0^2*α0*L/(4*π^2*σ12)
    Ppump = Pp*(1-exp(-α0*β*L*(1-N)))/(N0*π*r0^2*L)

    du[1] = 2*L/Tr*P*(rw*α0*(N*(ξ1 - ξ2) - 1) - αth) + Psp
    du[2] = -σ12*rw*P/(π*r0^2)*(N*ξ1 - 1) - N/τ + Ppump
end




function compute_basins_lasers(di::Dict)
    @unpack  res = di
    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e9)
    ds = CoupledODEs(pumped_laser!, rand(2), []; diffeq)
    yg = range(-1, 1; length = 10001)
    xg = range(-1, 30; length = 10001)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid)
    yr = range(0, 1, length = res)
    xr = range(0, 15, length = res)
    bsn = @showprogress [ mapper([x,y]) for x in xr, y in yr]
    grid = (xr,yr)
    return @strdict(bsn, grid, res)
end

res = 100
params = @strdict res
print_fig(params, "basins_lasers", compute_basins_lasers; force = true)

