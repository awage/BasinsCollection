using SparseArrays
using LsqFit
using Random 
function unwrap!(yphase)
    for i in 2:length(yphase)
      while yphase[i] - yphase[i-1] >= pi
        yphase[i] -= 2pi
      end
      while yphase[i] - yphase[i-1] <= -pi
        yphase[i] += 2pi
      end
    end
end

function proj_fun(y)
    @. model(x, p) = x*p[1]+p[2]
    p0 = [0.5, 0.5]
    xdata = 1:length(y)
    N = length(y)
    unwrap!(y)
    fit = curve_fit(model, xdata, y, p0)
    p = fit.param
    p[1] = p[1]/(2π/N)
    p[2] = 0.
    return p
end

function _proj_fun_f(y,t)
    p = proj_fun(Array(y[end]))
    return [p[1]] 
end


default_diffeq = (alg = Vern9(), reltol = 1e-9,  maxiters = 1e8)
mutable struct KuramotoParameters{M}
    N::Int
    Δ::M
    x::Vector{Float64}
    y::Vector{Float64}
end
function KuramotoParameters(; N)
    Δ = SparseMatrixCSC(Bidiagonal(-ones(N), ones(N-1), :U))
    Δ[N,1] = 1
    x = Δ * zeros(N)
    y = Δ * zeros(N)
    return KuramotoParameters(N, Δ, x, y)
end

# Optimized version of the ring Kuramoto network
# There is no allocation, the temporal variables are 
# stored in the parameters x and y 
function ring_kuramoto!(du, u, p, t)
    (; N, Δ, x, y) = p
    mul!(x, Δ, u) 
    x .= sin.(x)
    mul!(y, Δ', u)
    y .= sin.(y)
    @. du = x + y
    return nothing
end

function ring_kuramotos(; N=10, diffeq = default_diffeq)
    p = KuramotoParameters(; N)
    ds = CoupledODEs(ring_kuramoto!, zeros(N), p; diffeq)
    return ds
end

