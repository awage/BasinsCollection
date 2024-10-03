using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
# using OrdinaryDiffEq:Vern9
using OrdinaryDiffEq
using ProgressMeter
include(srcdir("print_fig.jl"))

# Typical bifurcation scenario in a three region identical New
# Economic Geography model
# Pasquale Commendatore , Ingrid Kubin, Iryna Sushko
# http://dx.doi.org/10.1016/j.matcom.2014.01.004
function f(Mr,Ms) 
    if Mr ≤ 0 
        return 0
    elseif Mr > 0 && Ms > 0 && Mr + Ms < 1 
        return Mr
    elseif Mr > 0 && Ms > 0 && Mr + Ms ≥ 1
        return Mr/(Mr + Ms)
    elseif Mr > 0 && Ms ≤ 0 && Mr + Ms < 1 
        return Mr/(1 - Ms)
    elseif Mr > 0 && Ms ≤ 0 && Mr + Ms ≥ 1
        return 1 
    end
end

function compute_aux_constants(λ1, λ2, λ3, μ, γ, σ, Φ)
    Δ1 = λ1 + Φ*(1-λ1); Δ2 = λ2 + Φ*(1-λ2); Δ3 = 1- (λ1 + λ2)*(1-Φ)
    C1 = (σ-μ)/3/(σ -μ*λ1*(1/Δ1 - Φ/Δ3))
    C2 = μ*Φ*λ1/(σ -μ*λ1*(1/Δ1 - Φ/Δ3))
    C3 = (σ-μ)/3/(σ -μ*λ2*(1/Δ2 - Φ/Δ3))
    C4 = μ*Φ*λ2/(σ -μ*λ2*(1/Δ2 - Φ/Δ3))
    # We solve the linear system for s1, s2, s3
    A = [1 -C2*(1/Δ2 - 1/Δ3) 0; 
         -C4*(1/Δ1 - 1/Δ3) 1 0;
         1 1 1]
    L = [C1+C2/Δ3; C3+C4/Δ3; 1]
    S = A\L 
    s1 = S[1]; s2 = S[2]; s3 = S[3]
        # s1 = ((σ-μ)/3 +μ*Φ*λ1*(s2/Δ2+(1-s2)/Δ3))/(σ -μ*λ1*(1/Δ1 - Φ/Δ3))
        # s2 = ((σ-μ)/3 +μ*Φ*λ2*(s1/Δ1+(1-s1)/Δ3))/(σ -μ*λ2*(1/Δ2 - Φ/Δ3))
        # s3 = 1-s1-s2
    D = λ1*Δ1^(μ/(σ-1))*(s1/Δ1 + Φ*(s2/Δ2 + s3/Δ3)) + λ2*Δ2^(μ/(σ-1))*(s2/Δ2 + Φ*(s1/Δ1 + s3/Δ3)) + λ3*Δ3^(μ/(σ-1))*(s3/Δ3 + Φ*(s1/Δ1 + s2/Δ2)) 
    K1 = Δ1^(μ/(σ-1))*(s1/Δ1 + Φ*(s2/Δ2 + s3/Δ3))/D
    K2 = Δ2^(μ/(σ-1))*(s2/Δ2 + Φ*(s1/Δ1 + s3/Δ3))/D
    M1 = λ1*(1+γ*(K1-1))
    M2 = λ2*(1+γ*(K2-1))
    return M1, M2
end

function economic_model!(du, u, p, t)
    μ = 0.4; γ = 5.; σ = 5.; Φ = p[1] 
    λ1, λ2 = u;  λ3 = 1 - (λ1 + λ2)
    M1, M2 = compute_aux_constants(λ1, λ2, λ3, μ, γ, σ, Φ)
    Z1 = f(M1,M2) 
    Z2 = f(M2,M1) 
    du[1] = Z1
    du[2] = Z2
    return nothing
end


function compute_basins_economic(di::Dict)
    @unpack  res,Φ = di
    ds = DeterministicIteratedMap(economic_model!, [0.5, 0.5], [Φ])
    yg = xg =  range(-1e-5, 1; length = 5000)
    grid = (xg, yg)
    mapper = AttractorsViaRecurrences(ds, grid; 
    consecutive_recurrences = 4000
    )
    xr = yr = range(1e-5, 1, length = res)
    bsn = zeros(res,res)
    @showprogress for (j,x) in enumerate(xr)
        for (k,y) in enumerate(yr)
            if x+y ≤ 1
                bsn[j,k] = mapper([x,y])
            else
              bsn[j,k] =  -1
            end
        end
    end
    grid = (xr,yr)
    att = extract_attractors(mapper)
    return @strdict(bsn, grid, res, att)
end


let res = 1500
    Φ = 0.085
    params = @strdict res Φ
    print_fig(params, "basins_economic", compute_basins_economic; force = false, xlab = L"\lambda_1", ylab = L"\lambda_2")
# att = get_att(params, "basins_economic", compute_basins_economic)
end

