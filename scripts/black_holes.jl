using DrWatson
@quickactivate # exports Attractors, GLMakie and other goodies in `src`
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ProgressMeter

# BASINS FOR MAJUMDAR-PAPAPETROU DOUBLE BLACK HOLE SYSTEM WITH SYMPLECTIC ALGORITHM*/
function double_black_hole_de!(dv, v, p, t)
    z1 = p[1]; z2 = p[2]; M = p[3]; pϕ = p[4] 
    z, pz, ρ, pρ = v
    t_ρ_z1 = (ρ^2+(z-z1)^2)
    t_ρ_z2 = (ρ^2+(z-z2)^2)
    t_M_z1_z2 = (1 + M/√t_ρ_z1 + M/√t_ρ_z2)
    ∂H∂z = pz/t_M_z1_z2^2
    ∂H∂pz = -2*t_M_z1_z2*((z-z1)*M/(t_ρ_z1^1.5) + (z-z2)*M/(t_ρ_z2^1.5))
    ∂H∂ρ = pρ/t_M_z1_z2^2
    ∂H∂pρ = -2*ρ*t_M_z1_z2*(M/(t_ρ_z1^1.5) + M/(t_ρ_z2^1.5)) + pϕ^2/ρ^3/t_M_z1_z2^2
    # ∂H∂ϕ = pϕ/ρ^2/(1+M/√t_ρ_z1+M/√t_ρ_z2)^2;
    dv[1] = ∂H∂z
    dv[2] = ∂H∂pz
    dv[3] = ∂H∂ρ
    dv[4] = ∂H∂pρ
end

function Hbh(p,q, prm) 
    z1 = prm[1]; z2 = prm[2]; M = prm[3]; pϕ = prm[4] 
    ρ, z = q;  pρ, pz = p
    U = 1 + M/sqrt(ρ^2 + (z- z1)^2) + M/sqrt(ρ^2 + (z - z2)^2)  
    h = ρ*U^2
    V = -0.5/ρ^2*(h-pϕ)*(h+pϕ)
    H = 0.5*(pρ^2 + pz^2) + V 
end

function salida_ρ_z(ρ, z, p)
    z1 = p[1]; z2 = p[2]; eps = 1/50
    sal =0 
    if √(ρ^2+z^2) > 100 # escaping photons
        sal=0;
    elseif abs(z-z1) + abs(ρ) < eps # upper shadow
        sal=1;
    elseif abs(z-z2) + abs(ρ) < eps # lower shadow
        sal=2;
    else 
        sal=-1;
    end    
return sal
end

function condition_ρ_z(ρ, z)
    z1 = 0.5 ; z2 = -0.5; eps = 1/50
    if √(ρ^2+z^2) > 100 # escaping photons
        return true
    elseif abs(z-z1) + abs(ρ) < eps # upper shadow
        return true
    elseif abs(z-z2) + abs(ρ) < eps # lower shadow
        return true
    else 
        return false
    end    
end


# Variables to store trajectories and energy evolution
# y0 = range(ρi, ρf, npuntos) 
# yp = range(zi, zf, npuntos) 
function get_basin_BH(ρ,z, prm)
    pϕ = prm[4]; M = prm[3]; z1 = prm[1]; z2 = prm[2]
    t_ρ_z1 = (ρ^2+(z-z1)^2)
    t_ρ_z2 = (ρ^2+(z-z2)^2)
    U = 1 + M/√t_ρ_z1 + M/√t_ρ_z2
    # @show U^4 - pϕ^2/ρ^2 
    if U^4 - pϕ^2/ρ^2 < 0
        return -1 
    else
        p = √(U^4-pϕ^2/ρ^2); # momento total p=√(p_ρ^2+p_z^2)
        pρ = p/√(1+(ρ-√3/2)^2/z^2); # tangential shooting from ρ=√(3)/2, z=0 (maximum of the potential) to determine p_ρ
        pz = √(p^2-pρ^2); # p_z from energy conservation
    end
    if  z > 0
        pρ = -pρ;
    end  
    if  ρ < √3/2 
        pz = -pz
    end
         
    vi = [z, pz, ρ, pρ];
    tspan=(0.0,50000.)
    # Define callback to halt solver
    condition(u,t,integrator) = condition_ρ_z(u[3], u[1])
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
    prob = ODEProblem(double_black_hole_de!, vi, tspan, prm, callback=cb)
    sol = solve(prob, Vern9(), reltol=1e-12, abstol=1e-12)
    u = sol[end]
    return salida_ρ_z(u[3], u[1], prm)
end


# get_basin_BH(1., 0.)



function compute_BH(di::Dict)
    @unpack Δpϕ, res = di
    φ = 0.5(1+√5)
    pϕ_crit = 0.5*5^(5/4)*φ^(3/2) # critical value 
    pϕ = pϕ_crit - Δpϕ
    M = 1; z1 = 0.5; z2 = -0.5        
    prm = [z1, z2, M, pϕ]
    zr = range(-0.7, 0.7; length = res)
    ρr = range(0., 2.; length = res)
    bsn = zeros(Int16, res, res)
    @showprogress for i in 1:res
        Threads.@threads for j in 1:res
            z = zr[j]; ρ = ρr[i];
            bsn[i,j] = get_basin_BH(ρ, z, prm) 
        end
    end
    grid = (ρr, zr)
    return @strdict(bsn, grid, Δpϕ, res)
end

Δpϕ = 0.03; #res = 500
params = @strdict Δpϕ res
print_fig(params, "black_holes", compute_BH; xlab = L"\rho", ylab = L"z", force = false)
