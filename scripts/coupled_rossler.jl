# using DrWatson
# @quickactivate "AdvancedBasins" # exports DynamicalSystems, GLMakie and other goodies in `src`
# using CairoMakie
# using LaTeXStrings
using DynamicalSystems
using OrdinaryDiffEq:Vern9
using Plots

# Reference: PHYSICAL REVIEW E 85, 035202(R) (2012)
#How to obtain extreme multistability in coupled dynamical systems
#C. R. Hens,1 R. Banerjee,1,2 U. Feudel,3 and S. K. Dana1
function coupled_rossler!(du, u, p, t)
    a = p[1]; b = p[2]; c = p[3];
    du[1] = - u[2] - u[3]
    du[2] = u[1] + a*u[5]
    du[3] = b - c*u[3] + u[4]*u[6]
    du[4] = u[1] - u[4] -u[2] - u[3]
    du[5] = u[4] + a*u[5]
    du[6] = b + u[6]*(u[4] -c)
end


a = 0.2; b = 0.2; c = 5.7
ds = ContinuousDynamicalSystem(coupled_rossler!, zeros(6), [a, b, c], (J,z0, p, n) -> nothing)
diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)


# _complete(y) = (length(y) == N) ? [Δϕ; Δω] : y; 
# _proj_state(y) = y[N+1:2*N]
# psys = projected_integrator(ds, _proj_state, _complete; diffeq)
yg = range(-25, 25; length = 10001)
grid = ntuple(x -> yg, dimension(ds))
mapper = AttractorsViaRecurrences(ds, grid; sparse = true, Δt = .1,   
    mx_chk_fnd_att = 100,
    mx_chk_loc_att = 100, safety_counter_max = Int(1e7), diffeq)
    #Ttr = 60.)

# Nt = 1000000
# l = [ mapper(5 .*(rand(6).-1)) for k in 1:Nt]
res = 100
y1 = range(-8, 8, length = res)
y2 = range(-2, 2, length = res)
ics = [ [-0.1, 0.01, 0.3, t1, t2, 2.] for t1 in y1, t2 in y2 ]
bas = [ mapper(u) for u in ics]
# for k = 1:Nt
#     @show l = mapper(5 .*(rand(6).-1)) 
# end
# u = trajectory(ds, 1000, [7.0047918247306695, -3.043913273644653, 0.21669156206552287, 7.0047918247306695, -5.056146159138251, 0.2166915620655228]; diffeq, Δt = 0.1)
