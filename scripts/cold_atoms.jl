using DrWatson
@quickactivate # exports Attractors, GLMakie and other goodies in `src`
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ProgressMeter

function cold_atoms_de!(dv,v,p,t)
    # α1 = p[1]; α2 = p[2]; β1 = p[3]; θ = p[4];   
    α1 = α2 = β1 = β2 = 1; θ = π/4;   
    x, px, y, py = v
    f(x,y) = (sin(θ)*x+cos(θ)*y)
    ∂H∂x = 2*α2*β2*sin(θ)*f(x,y)*exp(-β2*f(x,y)^2);
    ∂H∂y = 2*α1*β1*y*exp(-β1*y^2)+2*α2*β2*cos(θ)*f(x,y)*exp(-β2*f(x,y)^2);
    ∂H∂px = px;
    ∂H∂py = py;
    dv[1] = ∂H∂px;
    dv[2] = -∂H∂x;
    dv[3] = ∂H∂py;
    dv[4] = -∂H∂y;
end

function Hca(p,q, prm) 
    α1 = α2 = β1 = β2 = 1; θ = π/4;   
    x, y = q
    px, py = p
    H = 0.5*(px^2+py^2) - α1*exp(-β1*y^2) - α2*exp(-β2*(x*sin(θ) + y*cos(θ))^2)
end

function salida(sol)
    x = sol[3,end]
    y = sol[4,end]
    # @show (x,y, sol.t[end])
    if x > 20 && abs(y) < 4
        sal = 1;
    elseif x < -20 && abs(y) < 4
        sal = 3;
    elseif  y < -20
        sal = 4;
    elseif y > 20
        sal = 2
    else 
        sal = -1;
        # @show x,y
    end
    return sal
end


function get_cold_atoms(x, vx, y, vy)
    # α1 = p[1]; α2 = p[2]; β1 = p[3]; θ = p[4];   
    α1 = α2 = β1 = β2 = 1; θ = π/4;   
    T=0.5*(vx^2+vy^2); # Energia Cinetica Inicial
    V = -α1*exp(-β1*y^2)-α2*exp(-β2*(sin(θ)*x+cos(θ)*y)^2); # Energia Potencial Inicial
    if T+V > 0
        return -1
    end
    vi = [x, vx, y, vy];
    tspan=(0.0,50000.)
    # Define callback to halt solver
    condition(u,t,integrator) = (u[3]^2+u[4]^2) > 40^2 && t > 5500
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
    prob = ODEProblem(cold_atoms_de!, vi, tspan, callback=cb)
    sol = solve(prob, Vern9(), reltol=1e-12, abstol=1e-12)
    return salida(sol)
end

function compute_cold_atoms(di::Dict)
    @unpack vx, x, res = di
    yr = range(-3, 3, length = res)
    vyr = range(-1.5, 1.5, length = res)
    bsn = zeros(Int16, res, res)
    @showprogress for i in 1:length(yr)
        Threads.@threads for k in 1:length(vyr)
            bsn[i,k] = get_cold_atoms(x, vx, yr[i], vyr[k])
        end
    end
    grid = (yr, vyr)
    return @strdict(bsn, grid, x, vx, res)
end

x = -500; vx = 0.1; #res = 300
params = @strdict x vx res
cmap = ColorScheme([RGB(1,1,1), RGB(0.55,0.9,0.35), RGB(0.9,0.4,0.1),  RGB(0.50,0.24,1)] )
print_fig(params, "cold_atoms", compute_cold_atoms; ylab = L"v_y", xlab = L"y", cmap)
