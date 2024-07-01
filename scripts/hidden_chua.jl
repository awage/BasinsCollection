using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter
# using Plots

# Physics Letters A Volume 375, Issue 23, 6 June 2011, Pages 2230-2233
#Localization of hidden Chuaʼs attractors
# A.Leonov N.V. Kuznetsov  V.I.Vagaitsev
function chua!(du, u, p, t)
    α, β, γ, m0, m1 = p
    x,y,z = u
    ψ(x) = m1*x + 0.5*(m0-m1)*(abs(x+1)-abs(x-1))
    du[1] = α*(y - x - ψ(x))
    du[2] = x - y + z
    du[3] = -(β*y + γ*z)
end




function compute_basins_chua(di::Dict)
    @unpack α, β, γ, m0, m1, res = di

    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e9)
    ds = CoupledODEs(chua!, rand(3), [α, β, γ, m0, m1]; diffeq)
    yg = range(-5, 5; length = 10001)
    xg = range(-30, 30; length = 10001)
    zg = range(-30, 30; length = 10001)
    grid = (xg, yg, zg)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, Δt = 0.01,   
        mx_chk_hit_bas = 20000,
        mx_chk_fnd_att = 100,
        mx_chk_loc_att = 100, maximum_iterations = Int(1e8))
    yr = range(-5, 5, length = res)
    xr = range(-7, 7, length = res)
    bsn = @showprogress [ mapper([x,y,0.]) for x in xr, y in yr]
    grid = (xr,yr)
    return @strdict(bsn, grid, α, β, γ, m0, m1, res)
end



# α = 9.3515908493; β = 14.7903198054; γ = 0.0160739649; m0 = -1.1384111956; m1 = -0.7224511209;
# α = 8.5; β = 14.28; γ = 0; m0 = -8/7; m1 = -5/7
α = 8.4562218418; β = 12.0732335925; γ = 0.0051631393; m0 = -0.1767573476; m1 = -1.1467573476; # res = 500
params = @strdict α β γ m0 m1 res
print_fig(params, "hidden_chua", compute_basins_chua)

