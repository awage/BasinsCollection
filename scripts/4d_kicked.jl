using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using Attractors
using Attractors
using StaticArrays
include("../src/print_fig.jl")


# REFERENCE: Feudel, U., Grebogi, C., Poon, L., & Yorke, J. A. (1998). Dynamical properties of a simple mechanical system with a large number of coexisting periodic attractors. Chaos, Solitons & Fractals, 9(1-2), 171-180.
function kicked_feudel_4d!(dz, z, p, n)
    θ1 = z[1]; θ2 = z[2]; θd1 = z[3]; θd2 = z[4];
    M = p[1]; L = p[2]; ρ = p[3]
    m1 = m2 = 1.; l1 = 1/√2; l2 = 1.; 
    I1 = (m1 + m2)*l1^2; I2 = m2*l2^2
    dz[1] = rem2pi(θ1 + M[1,1]*θd1 + M[1,2]*θd2, RoundNearest) 
    dz[2] = rem2pi(θ2 + M[2,1]*θd1 + M[2,2]*θd2, RoundNearest) 
    dz[3] = ρ*l1/I1*sin(dz[1]) + L[1,1]*θd1 + L[1,2]*θd2
    dz[4] = ρ*l2/I2*sin(dz[2]) + L[2,1]*θd1 + L[2,2]*θd2
    return 
end



function compute_kicked_rotor_4d(di)
    @unpack M, L , ρ, res = di
    u0 = [0.; 0.; 0.; 0.]
    p = [M, L, ρ]
    df = DiscreteDynamicalSystem(kicked_feudel_4d!, u0, p) 
    θ1 = range(-2π, 2π, length = 10001)
    θ2 = range(-2π, 2π, length = 10001)
    grid_rec = (θ1, θ2)
    projection = [1,2]; complete_state = [9. , 2.5]
    pinteg = projected_integrator(df, projection, complete_state)
    mapper = AttractorsViaRecurrences(pinteg, grid_rec,
            # mx_chk_lost = 10, 
            mx_chk_fnd_att = 300, 
            mx_chk_loc_att = 300, 
            # mx_chk_att = 2,
             sparse = true
    )
    θ1 = range(-π, π, length = res)
    θ2 = range(-π, π, length = res)
    grid = (θ1, θ2)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, M, L, ρ, res)
end


 

ν = 0.2; ρ = 6.5; T = 1
ν1 = ν2 = ν
Δ = √(ν1^2 + 4*ν2^2) 
a = 0.5*(1 + ν1/Δ); d = 0.5*(1 - ν1/Δ); b = -ν2/Δ;
λ1 = -0.5*(ν1 + 2*ν2 + Δ) ; λ2 = -0.5*(ν1 + 2*ν2 - Δ) ; 
W1 = [a b; b d]; W2 = [d -b; -b a];
L = W1*exp(λ1*T) + W2*exp(λ2*T)  
M = W1*(exp(λ1*T)-1)/λ1 + W2*(exp(λ2*T) -1)/λ2
res = 400
params = @strdict res  M L  ρ
print_fig(params, "kicked_rotor_4d", compute_kicked_rotor_4d; ylab = L"\dot{\theta}", xlab = L"\theta")
