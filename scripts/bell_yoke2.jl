using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes
using ProgressMeter
include(srcdir("print_fig.jl"))

# Brzeski, P., Kapitaniak, T., & Perlikowski, P. (2015). Analysis of transitions between different ringing schemes of the church bell. International Journal of Impact Engineering, 85, 57-66.
function bell_yoke!(du, u, p, t)
    ϕ1, ϕ2, dϕ1, dϕ2 = u
    Tmax, lr = p
    M = 2633; m = 57.4; 
    Bb = 1375; Bc = 45.15 
    L = 0.236; l = 0.739; lc = -0.1 
    α = 0.5349; Db = 26.68; Dc = 4.539 
    g = 9.8; Lr = L - lr; lcr = lc -lr
    Bbr = (Bb - M*L^2) + M*Lr^2
    A = 15.; ω = 7.5 
    Mt(ϕ,dϕ) = abs(ϕ) ≤ π/A ? Tmax*sign(dϕ)*cos(ω*ϕ) : 0.

    du[1] = dϕ1; du[2] = dϕ2; 

    # AA = [(Bbr + m*lcr^2)  m*lcr*l*cos(ϕ2 - ϕ1);
    #       m*lcr*l*cos(ϕ2-ϕ1)  Bc]

    γ = m*lcr*l*cos(ϕ2 - ϕ1)
    AA_inv = 1/((Bbr + m*lcr^2)*Bc - γ^2 )*[ Bc -γ; -γ (Bbr + m*lcr^2)]
    b1 = m*lcr*l*dϕ2^2*sin(ϕ2-ϕ1) - (M*Lr + m*lcr)*g*sin(ϕ1) - Db*dϕ1 + Dc*(dϕ2-dϕ1) + Mt(ϕ1,dϕ1)
    b2 = - m*lcr*dϕ1^2*sin(ϕ2 - ϕ1) - m*g*sin(ϕ2) - Dc*(dϕ2 - dϕ1) 
    # du[3:4] = AA\[b1;b2]
    du[3:4] = AA_inv*[b1;b2]
end


# We have to define a callback to wrap the phase in [-π,π]
function affect!(integrator)
    uu = integrator.u
    ϕ1, ϕ2, dϕ1_BI, dϕ2_BI = uu
    M = 2633; m = 57.4; 
    Bb = 1375; Bc = 45.15 
    L = 0.236; l = 0.739; lc = -0.1 
    k =  0.05
    Lr = L - lr; lcr = lc -lr
    Bbr = (Bb - M*L^2) + M*Lr^2

    # Solve linear system
    a1 = (Bbr + m*lcr^2 + m*lcr*l*cos(ϕ2 - ϕ1))
    a2 = (Bc + m*lcr*l*cos(ϕ2 - ϕ1))
    BB = [a1  a2; 
          1 -1 ]
    bb1 = a1*dϕ1_BI + a2*dϕ2_BI
    bb2 = sqrt(k)*(dϕ2_BI - dϕ1_BI)
    dϕ_AI = BB\[bb1; bb2]
    dϕ1_AI, dϕ2_AI = dϕ_AI

    set_state!(integrator, SVector(uu[1], uu[2], dϕ1_AI, dϕ2_AI))
    u_modified!(integrator, true)
end

    # α = 0.5349 
    condition(u,t,integrator) = abs(u[1] - u[2]) - 0.5349

function compute_bell_yoke(di::Dict)
    @unpack Tmax, lr, res = di
    cb = ContinuousCallback(condition,affect!; abstol = 1e-18, save_positions=(true,true), interp_points = 30)
    diffeq = (reltol = 0, abstol = 1e-17,  alg = Vern9(), callback = cb)
    df = CoupledODEs(bell_yoke!, rand(4), [Tmax, lr]; diffeq)
    # psys = ProjectedDynamicalSystem(df, [1,2], [0.0, 0.0])
    xg = range(-0.7, 0.7; length = 101)
    yg = zg = tg = range(-2, 2; length = 101)
    grid = (xg, yg, zg, tg)
    mapper = AttractorsViaRecurrences(df, grid; 
            Ttr = 200,
            # maximum_iterations = 1e6, 
            consecutive_attractor_steps = 2)
            # consecutive_basin_steps = 200,
            # consecutive_recurrences = 200)
    xg = range(-0.5,0.5,length = res)
    yg = range(-0.5,0.5,length = res)
    grid = (xg, yg)
    bsn = zeros(Int, res, res)
@showprogress    for (i,x) in enumerate(xg), (j,y) in enumerate(yg)
        if abs(x - y) < 0.5349
            bsn[i,j] = mapper([x, y, 0., 0.])
        end
    end
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid, res)
end


res = 700; Tmax = 150; lr = -0.03
params = @strdict res Tmax lr
cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
print_fig(params, "bell_yoke", compute_bell_yoke; ylab= L"\phi_1", xlab= L"\phi_2", force = false, cmap)
att = get_att(params, "bell_yoke", compute_bell_yoke)


# Tmax = 125; lr = 0.20
# Tmax = 175; lr = 0.16
# Tmax = 325; lr = -1.21
# Tmax = 150; lr = -0.03
# α = 0.5349 
# Δ = 0.0001
# cond(u,t,integrator) = abs(u[1] -  u[2]) > α + Δ
# function aff!(integrator)
#     uu = integrator.u
#     # println("discrete callback")
#     if uu[1] - uu[2] > α 
#         # set_state!(integrator, SVector(uu[1], uu[1] + α - Δ/2, uu[3], uu[4]))
#         set_state!(integrator, SVector(uu[2] + α - Δ/2 ,uu[2], uu[3], uu[4]))
#         u_modified!(integrator, true)
#     elseif uu[2] - uu[1] > α
#         # set_state!(integrator, SVector(uu[1],uu[1] + α - Δ/2, uu[3], uu[4]))
#         set_state!(integrator, SVector(uu[2] - α + Δ/2, uu[2], uu[3], uu[4]))
#         u_modified!(integrator, true)
#     end
# end

# cbd = DiscreteCallback(cond, aff!)
# cb = ContinuousCallback(condition,affect!; abstol = 1e-18, save_positions=(true,true), interp_points = 30)
# diffeq = (reltol = 0, abstol = 1e-17,  alg = Vern9(), callback = CallbackSet(cb,cbd))
# df = CoupledODEs(bell_yoke!, rand(4), [Tmax, lr]; diffeq)
# y,t = trajectory(df, 2000, rand(4)*0.05)
# # lines(y[:,1],y[:,2])
# # lines(y[19000:20000,1],y[19000:20000,2])
# p = Figure(); ax = Axis(p[1,1])
# lines!(ax, [-0.75; -0.2],  [-0.75+α; -0.2+α])
# lines!(ax,y[18000:18200,1],y[18000:18200,2])
# p
