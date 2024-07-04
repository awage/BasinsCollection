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
    # α, β, γ, m0, m1 = p
    # L = 70; Tr = 8.7; r0 = 1.4e-4; σ12 = 3e-21; 
    # τ = 1e-2; ξ1 = 2; ξ2= 0.4; N0 = 1.3e20; R = 0.8
    # w0 = 3.5e-4; γ0 = 0.038; λg = 1.56e-4; β = 0.5
    # α0 = N0*σ12 
    # αth = γ0 + 1/(2*L)*log(1/R)
    # # rw = 1 - exp(-2*(r0/w0)^2)
    # rw = 0.308
    # P,N = u
    md = 0.5; fd = 35e3
    # md = 0.8; fd = 70.2e3
    Pp = 7.4e19*(1+md*sin(2*π*fd*t/2.8e8))
    # Psp = N*1e-3/(τ*Tr)*(λg/w0)^2*r0^2*α0*L/(4*π^2*σ12)
    # Ppump = Pp*(1-exp(-α0*β*L*(1-N)))/(N0*π*r0^2*L)

    # du[1] = 2*L/Tr*P*(rw*α0*(N*(ξ1 - ξ2) - 1) - αth) + Psp
    # du[2] = -σ12*rw*P/(π*r0^2)*(N*ξ1 - 1) - N/τ + Ppump
    
    a1 = 2.4; a2 = 6.9e-13; a3 = 5.1e-13; b1 = 3.5e-7; 
    b2 = 2.6e-7; b3 = 0.5; P0 = 2e-23*Pp
    du[1] = u[1]*u[2] - a1*u[1] + a2*u[2] + a3
    du[2] = -u[1]*u[2] - b1*u[2] - b2 + P0*(1-b3*exp(u[2])) 
end




function compute_basins_lasers(di::Dict)
    @unpack  res = di

    diffeq = (alg = Vern9(), reltol = 1e-11, maxiters = 1e9)
    ds = CoupledODEs(pumped_laser!, rand(2), []; diffeq)
    xg = range(0, 20; length = 20000)
    yg = range(0, 1; length = 20000)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid;    
        consecutive_recurrences= 10000 )
    yr = range(0.01, 10, length = res)
    xr = range(0.01, 1, length = res)
    bsn = @showprogress [ mapper([x,y]) for x in xr, y in yr]
    grid = (xr,yr)
    return @strdict(bsn, grid, res)
end


res = 100
params = @strdict res
print_fig(params, "basins_lasers", compute_basins_lasers; force = true)

