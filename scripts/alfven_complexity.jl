using DrWatson
@quickactivate # exports Attractors, GLMakie and other goodies in `src`
using Attractors
using OrdinaryDiffEq:Vern9
using CairoMakie
using LaTeXStrings

# ALFVÉN BOUNDARY CRISIS
# A. C.-L. CHIAN, F. A. BOROTTO and E. L. REMPEL
# https://doi.org/10.1142/S0218127402005303
@inline @inbounds function alfven(u, p, t)
    a = 0.1; ω = -1; λ = 1/4; ν = p[1]
    ∂Hy(y,z) = (z^2 + y^2 -1)*y - λ*(y-1)
    ∂Hz(y,z) = (z^2 + y^2 -1)*z - λ*z
    du1 = ∂Hz(u[1],u[2]) + a*cos(ω*t)
    du2 = -∂Hy(u[1],u[2]) + a*sin(ω*t)
    du1, du2 = inv([1 -ν; ν 1])*[du1; du2]
    return SVector{2}(du1, du2)
end

function compute_alfven(di::Dict)
    @unpack  res,ν = di
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    ds = CoupledODEs(alfven, rand(2), [ν]; diffeq)
    xg = yg = range(-3,3,length = 50000)
    ω = 1
    smap = StroboscopicMap(ds, 2*pi/ω)
    mapper = AttractorsViaRecurrences(smap, (xg, yg); 
       consecutive_recurrences = 1000 )
    xg = yg = range(-2.2,2.2,length = res)
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid,  res)
end

let res = 1200
ν = 0.01747; 
params = @strdict ν res
print_fig(params, "alfven", compute_alfven; ylab = L"$\dot{x}$", xlab = L"x", force = false) 
end

