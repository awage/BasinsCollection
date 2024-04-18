using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using Attractors
using OrdinaryDiffEq:Vern9
using ProgressMeter
include("../src/print_fig.jl")
#
# CHAOS 28, 113110 (2018)
# Multistability in the cyclic competition system
# Junpyo Park,1 Younghae Do,2 and Bongsoo Jang
# https://doi.org/10.1063/1.5045366
function cyc_comp!(du, u, p, t)
    r = p; μ = σ = 1
    a,b,c = u; ρ = a + b + c
    du[1] = a*(μ*(1-ρ) - σ*c - r*a^2*(1-a)/2)
    du[2] = b*(μ*(1-ρ) - σ*a - r*b^2*(1-b)/2)
    du[3] = c*(μ*(1-ρ) - σ*b - r*c^2*(1-c)/2)
end




function compute_cyc_comp(di::Dict)
    @unpack r, res = di

    ds = ContinuousDynamicalSystem(cyc_comp!, rand(3), r, (J,z0, p, n) -> nothing)
    diffeq = (alg = Vern9(), reltol = 1e-10, maxiters = 1e9)
    yg = xg = zg = range(0., 4; length = 10001)
    
    mapper = AttractorsViaRecurrences(ds, (xg, yg, zg); sparse = true, Δt = .1,
        mx_chk_att = 2,
        mx_chk_fnd_att = 1000,
        mx_chk_loc_att = 1000, mx_chk_safety = Int(1e8), diffeq)
    # We force the mapper to find the attractors first. For some reasons it cannot find them all normally. 
    mapper([sqrt(2/r), 0, 0]); mapper([0, sqrt(2/r), 0]); mapper([0, 0, sqrt(2/r)])

    # Projection on a the plane x + y + z = \sqrt(2/r). u1 and u2 are orthogonal vectors on this plane. The function ics produces initial conditions on this plane.    
    u1 = [-1, 1, 0]/sqrt(2);  u2 = [-1, -1, 2]/sqrt(6); 
    zr = range(0., sqrt(3/r), length = res)
    yr = range(0., sqrt(4/r), length = res)
    function ics(y,z)         
        u = y*u1 + z*u2 + [sqrt(2/r), 0, 0]; #@show u
        if abs(sum(u) - sqrt(2/r)) > 0.01 || u[1] < 0. || u[2] < 0. 
            return -1
        else
            return mapper(u)
        end
    end    

    bsn = @showprogress [ics(y,z) for y in yr, z in zr]
    att = mapper.bsn_nfo.attractors
    grid = (yr,zr)
    return @strdict(bsn, grid, att, res)
end


r = 3.35; res = 1500;
params = @strdict r res
print_fig(params , "cyc_comp", compute_cyc_comp; xlab = "", ylab = "") 

