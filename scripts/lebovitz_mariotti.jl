using DrWatson
@quickactivate
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9
using CairoMakie


# Lebovitz, N., & Mariotti, G. (2013). Edges in models of shear flow. Journal of Fluid Mechanics, 721, 386-402. doi:10.1017/jfm.2013.38
mutable struct LMD6Parameters{M}
    k::M
    σ_0::Float64 
    σ_1::Float64
    σ_2::Float64
    σ_3::Float64
    σ_4::Float64
    σ_5::Float64
    σ_6::Float64
    Re::Float64
end
function LMD6Parameters(; Re = 307.)
    α = 1.1; β = π/2; γ = 5/3; 
    c_3 = sqrt(4*β^2 + γ^2); 
    c_5 = sqrt(γ^2 + α^2); 
    c_6 = sqrt(α^2*β^2 + γ^4 + 2*γ^2*α^2 + α^4 + 3/4*β^2*γ^2); 
    k = [   β^2; 
            5*β^2+ γ^2;
            c_3^2; 
            α^2 + 4*β^2; 
            α^2 + β^2 + γ^2; 
            (α^2 + β^2) + (γ^2*(4*c_5^4 + β^2*(4*α^2+γ^2)))/c_6^2]
    σ_0 = β*γ/c_3; 
    σ_1 = γ^2/c_5; 
    σ_2 = α^2/c_5; 
    σ_3 = γ*α*β/(2*c_6); 
    σ_4 = β^2*(4*α^2+5*γ^2)*α/(2*c_3*c_5*c_6)
    σ_5 = (β^2 - α^2 - γ^2)*γ^2*α/(2*c_3*c_5*c_6)
    σ_6 = γ^2*β^2*α/(4*c_3*c_5*c_6)

    return LMD6Parameters(k, σ_0, σ_1, σ_2, σ_3, σ_4, σ_5, σ_6, Re)
end

function LMD6!(du, u, p, t)
    (; k, σ_0, σ_1, σ_2, σ_3, σ_4, σ_5, σ_6, Re) = p
    du[1] = -u[1]*k[1]/Re -σ_0*u[2]*u[3]; 
    du[2] = -u[2]*k[2]/Re + σ_0*u[3] + σ_0*u[1]*u[3] - σ_1*u[4]*u[5]; 
    du[3] = -u[3]*k[3]/Re - (σ_4 + σ_5)*u[5]*u[6]; 
    du[4] = -u[4]*k[4]/Re - σ_3*u[6] + σ_2*u[2]*u[5] - σ_3*u[1]*u[6]; 
    du[5] = -u[5]*k[5]/Re + (σ_1 - σ_2)*u[2]*u[4] + (σ_4 - σ_6)*u[3]*u[6]; 
    du[6] = -u[6]*k[6]/Re + σ_3*u[4] + (σ_5 + σ_6)*u[3]*u[5] + σ_3*u[1]*u[4] 
end

function compute_LM(di::Dict)
    @unpack res, Re = di
    p = LMD6Parameters(; Re = Re)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    ds = CoupledODEs(LMD6!, zeros(6), p, diffeq)
    yg = range(-5, 5; length = 10001)
    grid = ntuple(x -> yg, dimension(ds))
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, Δt = .1,   
        mx_chk_fnd_att = 300,
        mx_chk_loc_att = 100, maximum_iterations = Int(1e7), diffeq)
    u0(x,y) = [x, -0.0511, -0.0391, 0.0016, y, 0.126]
    y1r = range(-1, 1, length = res)
    y2r = range(-1, 1, length = res)
    ics = [ u0(y1,y2) for y1 in y1r, y2 in y2r]
    bsn = [ mapper(u) for u in ics]
    grid = (y1r,y2r)
    return @strdict(bsn, grid, Re, res)
end


# res = 450; 
Re = 307.
params = @strdict res Re
print_fig(params, "lebovitz", compute_LM)
