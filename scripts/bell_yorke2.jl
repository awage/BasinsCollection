using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes
using ProgressMeter
include(srcdir("print_fig.jl"))

# ( B_b + m*lc^2 )*ddϕ1 + m*lc*ddϕ2 cos(ϕ2 - ϕ1)
# B_c*ddϕ2 + m*lc*ddϕ1*cos(ϕ2 - ϕ1) 

# AA = [( B_b + m*lc^2 )  m*lc*cos(ϕ2 - ϕ1);
# m*lc*cos(ϕ2 - ϕ1) B_c  ]

#  b1 = m*lc*dϕ2^2*sin(ϕ2 - ϕ1) - (ML + mlc)*g*sin(ϕ1) - Db*dϕ1 + Dc*(dϕ2 - dϕ1) + T(ϕ1)
# b2 = - m*lc*dϕ1^2*sin(ϕ2 - ϕ1) - m*g*sin(ϕ2) - Dc*(dϕ2 - dϕ1) 


# # (dϕ2_AI - dϕ1_AI) = sqrt(k)*(dϕ2_BI - dϕ1_BI)
# k = 0.05
# BB = [( Bb + ml^2_c + m*lc*l*cos(ϕ2 - ϕ1))  (Bc + m*lc*l*cos(ϕ2 - ϕ1)) ; 
# 1 -1 ]

# bb1 = (Bb + m*lc^2 + m*lc*l*cos(ϕ2 - ϕ1))*dϕ1_BI + (Bc + m*lc*l*cos(ϕ2 - ϕ1))*dϕ2_BI
# bb2 = sqrt(k)*(dϕ2_BI - dϕ1_BI)


# M = 2633 
# m = 57.4 
# Bb = 1375 
# Bc = 45.15 
# L = 0.236 
# l = 0.739 
# lc = -0.1 
# α = 0.5349 
# Db = 26.68 
# Dc = 4.539 
# Tmax = 229.6
# # lc =  0.25


function bell_yoke!(du, u, p, t)
    M = 2633 
    m = 57.4 
    Bb = 1375 
    Bc = 45.15 
    L = 0.236 
    l = 0.739 
    lc = -0.1 
    α = 0.5349 
    Db = 26.68 
    Dc = 4.539 
    g = 9.8
    Tmax = 150
    lr = -0.03
    Lr = L - lr
    lcr = lc -lr
    Bbr = (Bb - M*L^2) + M*Lr^2

    A = 15.; ω = 7.5 
    Mt(ϕ,dϕ) = abs(ϕ) ≤ π/A ? Tmax*sign(dϕ)*cos(ω*ϕ) : 0.

    ϕ1, ϕ2, dϕ1, dϕ2 = u
    du[1] = dϕ1; du[2] = dϕ2; 

    AA = [(Bbr + m*lcr^2)  m*lcr*l*cos(ϕ2 - ϕ1);
          m*lcr*l*cos(ϕ2-ϕ1)  Bc]

    b1 = m*lcr*l*dϕ2^2*sin(ϕ2-ϕ1) - (M*Lr + m*lcr)*g*sin(ϕ1) - Db*dϕ1 + Dc*(dϕ2-dϕ1) + Mt(ϕ1,dϕ1)
    b2 = - m*lcr*dϕ1^2*sin(ϕ2 - ϕ1) - m*g*sin(ϕ2) - Dc*(dϕ2 - dϕ1) 
    du[3:4] = AA\[b1;b2]
end


# We have to define a callback to wrap the phase in [-π,π]
function affect!(integrator)
    uu = integrator.u
    ϕ1, ϕ2, dϕ1_BI, dϕ2_BI = uu
    M = 2633 
    m = 57.4 
    Bb = 1375 
    Bc = 45.15 
    L = 0.236 
    l = 0.739 
    lc = -0.1 
    α = 0.5349 
    Db = 26.68 
    Dc = 4.539 
    lr = -0.03
    Lr = L - lr
    lcr = lc -lr
    Bbr = (Bb - M*L^2) + M*Lr^2
    a1 = (Bbr + m*lcr^2 + m*lcr*l*cos(ϕ2 - ϕ1))
    a2 = (Bc + m*lcr*l*cos(ϕ2 - ϕ1))
    k =  0.05
    BB = [a1  a2; 
          1 -1 ]
    bb1 = a1*dϕ1_BI + a2*dϕ2_BI
    bb2 = sqrt(k)*(dϕ2_BI - dϕ1_BI)
    dϕ_AI = BB\[bb1; bb2]
    dϕ1_AI, dϕ2_AI = dϕ_AI
    # @show dϕ1_BI, dϕ2_BI, dϕ_AI
    # @show (dϕ1_BI - dϕ2_BI)^2
    # @show (dϕ_AI[1] - dϕ_AI[2])^2
    # @show a1*dϕ1_BI + a2*dϕ2_BI
    # @show a1*dϕ1_AI + a2*dϕ2_AI
    set_state!(integrator, SVector(uu[1], uu[2], dϕ1_AI, dϕ2_AI))
    # set_state!(integrator, SVector(uu[1], uu[2], -sign(dϕ1_BI)*abs(dϕ_AI[1]), -sign(dϕ2_BI)*abs(dϕ_AI[2])))
    # set_state!(integrator, SVector(uu[1], uu[2], uu[3], -uu[4]))
    u_modified!(integrator, true)
end

    # α = 0.5349 
    condition(u,t,integrator) = abs(u[1] - u[2]) - 0.5349

function compute_bell_yoke(di::Dict)
    @unpack res = di
    cb = ContinuousCallback(condition,affect!)
    diffeq = (reltol = 0, abstol = 1e-14,  alg = Vern9(), callback = cb)
    df = CoupledODEs(bell_yoke!, rand(4); diffeq)
    # comp_state(y) = (length(y) == 2) ? rand(8) : y; 
    # psys = ProjectedDynamicalSystem(df, [1,2], comp_state)
    yg = range(-2, 2; length = 1001)
    grid = ntuple(x -> yg, 4)
    mapper = AttractorsViaRecurrences(df, grid; maximum_iterations = 1e6, Ttr = 100)
    xg = range(-0.5,0.5,length = res)
    yg = range(-0.5,0.5,length = res)
    grid = (xg, yg)
    # bsn = @showprogress [mapper(rand(4)*0.2) for x in xg, y in yg]  
    bsn = @showprogress [mapper([x, y, 0., 0.]) for x in xg, y in yg]  
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid, res)
end


res = 200
params = @strdict res
print_fig(params, "bell_yoke", compute_bell_yoke; ylab= L"\dot{\theta}", xlab= L"\theta", force = true)
att = get_att(params, "bell_yoke", compute_bell_yoke)

# cb = DiscreteCallback(condition,affect!)
# cb = ContinuousCallback(condition,affect!)
# diffeq = (reltol = 0, abstol = 1e-16,  alg = Vern9(), callback = cb)
# # diffeq = (reltol = 1e-6,  alg = Vern9(), maxiter = 1e6)
# df = CoupledODEs(bell_yoke!, rand(4); diffeq)
# # diffeq = (reltol = 1e-8,  alg = Vern9(),)
# # df = CoupledODEs(bell_yoke!, 5*rand(4); diffeq)
# y,t = trajectory(df, 250, rand(4)*0.05)
# # plot(t,y[:,1])
# lines(y[:,1],y[:,2])
