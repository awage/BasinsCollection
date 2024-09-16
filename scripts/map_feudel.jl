using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9



# Basin bifurcation in quasiperiodically forced systems Ulrike Feudel, Annette Witt, Ying-Cheng Lai, and Celso Grebogi PRE 28, 1998
# https://doi.org/10.1103/PhysRevE.58.3060
function chaotic_map(dz, z, p, n)
    xn = z[1]
    θ = z[2]
    a = p[1]
    ω = (sqrt(5.0) - 1.0) / 2.0
    r = p[2]
    f(x) = r * x * (1.0 - x)
    Mn(n) = reduce(∘, fill(f, n))
    M = Mn(3)
    dz[1] = M(xn) + a * cos(2 * π * θ)
    dz[2] = mod(θ + ω, 1.0)
    return
end


function compute_feudel(di::Dict)
    @unpack a, r, res = di
    ds = DeterministicIteratedMap(chaotic_map, [1.0, 0.0], [a, r])
    θ = range(0.0, 1.0, length = 2500)
    xg = range(0.0, 1.0, length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,θ))
    θ = range(0.0, 1.0, length = res)
    xg = range(0.0, 1.0, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,θ); show_progress = true)
    grid = (xg, θ)
    return @strdict(bsn, att, grid, b, res)
end


r = 3.833
a = 0.0015
# res = 400
params = @strdict res a r
print_fig(params, "feudel", compute_feudel; xlab = L"x", ylab = L"\theta", force = false)
