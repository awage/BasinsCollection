using DrWatson
@quickactivate
using LaTeXStrings
using CairoMakie
using ProgressMeter
using Random
include(srcdir("print_fig.jl"))

# Arithmetic of the sync basin for pulse-coupled oscillators
# K. P. O'Keeffe
# arXiv:2609.01668 (2026)
#
# N identical oscillators on a complete graph with voltage xᵢ ∈ [0,1]:
#       ẋᵢ = 1 - γ xᵢ
# When a cluster of size k reaches xᵢ = 1 it fires, resets to 0 and
# delivers a pulse Δ = k/N to all the others. Any oscillator pushed to
# xⱼ + Δ ≥ 1 is absorbed into the firing cluster (no cascade: only the
# firer's pulse counts). The number of clusters never increases, so the
# system ends in full synchrony or in a frozen state of coexisting clusters.
# At γ = 0 the surviving states are equal-size clusters whose size divides N.
#
# The dynamics is event driven: between firings the flow is solved exactly
# with w = 1 - γx, w(t) = w(0) exp(-γt), and we jump from firing to firing.

# Time for an oscillator at voltage x to reach the threshold x = 1.
time_to_fire(x, γ) = γ == 0 ? 1 - x : log((1 - γ*x)/(1 - γ))/γ
# Free evolution of the voltage x during a time t.
evolve(x, t, γ) = γ == 0 ? x + t : (1 - (1 - γ*x)*exp(-γ*t))/γ

# Run the event-driven dynamics from the voltages x0. Returns the sorted
# vector of final cluster sizes ([N] means synchrony). The run stops when a
# single cluster is left or when no absorption happened during `nrounds`
# full firing rounds (the cluster state is then considered frozen).
function pco_cluster_sizes(x0, γ; nrounds = 200, tol = 1e-12)
    N = length(x0)
    x = collect(float.(x0))
    k = ones(Int, N)            # cluster sizes
    quiet = 0                   # firings since the last absorption
    while length(x) > 1 && quiet < nrounds*length(x)
        f = argmax(x)
        t = time_to_fire(x[f], γ)
        @. x = evolve(x, t, γ)
        Δ = k[f]/N
        absorbed = false
        keep = trues(length(x))
        for j in eachindex(x)
            j == f && continue
            x[j] += Δ
            if x[j] ≥ 1 - tol
                k[f] += k[j]; keep[j] = false; absorbed = true
            end
        end
        x[f] = 0.
        x = x[keep]; k = k[keep]
        quiet = absorbed ? 0 : quiet + 1
    end
    return sort(k)
end

# Label of the final state: the number of clusters (1 = synchrony).
pco_outcome(x0, γ; kwargs...) = length(pco_cluster_sizes(x0, γ; kwargs...))


# Two dimensional slice of the N dimensional basin: the voltages x₁ and x₂
# span [0,1]² while the other N-2 voltages are kept fixed at random values.
function compute_pco_basins(di::Dict)
    @unpack N, γ, res = di
    rng = Xoshiro(1234)
    xr = rand(rng, N - 2)
    x1 = range(0, 1, length = res)
    x2 = range(0, 1, length = res)
    bsn = zeros(Int, res, res)
    p = Progress(res)
    Threads.@threads for i in 1:res
        for j in 1:res
            bsn[i,j] = pco_outcome([x1[i]; x2[j]; xr], γ)
        end
        next!(p)
    end
    grid = (x1, x2)
    return @strdict(bsn, grid, res, xr)
end


# Convex charging: unequal clusters (e.g. (4,3) for N = 7) can survive.
let res = 1200, N = 117, γ = -1.0
    params = @strdict res N γ
    print_fig(params, "pco_arithmetic", compute_pco_basins; force = false, xlab = L"x_1", ylab = L"x_2")
end
