using DrWatson
@quickactivate
using LaTeXStrings
using CairoMakie
using ProgressMeter
using LinearAlgebra
using Sundials
include(srcdir("print_fig.jl"))

# The basins of attraction of soft sphere packings are not fractal
# P. Suryadevara, M. Casiulis, S. Martiniani
# arXiv:2409.12113 (2026), Fig. 1 (CVODE panel)
# Authors' code: https://github.com/martiniani-lab/basinerror
#
# N polydisperse disks in a periodic square box of side L interacting with
# the Hertzian repulsion
#       V(r) = (1 - r/σ)^(5/2)   for r < σ = Rᵢ + Rⱼ
# Half the radii are drawn from N(1.0, 0.05), the other half from N(1.4, 0.07).
# L is fixed by the packing fraction ϕ = π ∑Rᵢ²/L². The state X ∈ R^(2N) follows
# the steepest descent flow
#       Ẋ = -∇E(X)
# integrated with CVODE (BDF) until ‖∇E‖ < gtol. The basins live on the 2d
# plane Xc + p n̂₁ + q n̂₂ with (p, q) ∈ [-6, 6] × [-3.375, 3.375].
# Two minima are the same if they share the rattlers and, after superposing the
# first non-rattler particle, all other non-rattler particles are closer than dtol.

# Configuration of the authors for N = 16, ϕ = 0.9 and seed 0, generated with
# numpy as in basinerror/utils/params.py (radii and X₀ from default_rng(0)) and
# basinerror/ensemble/random_plane.py (n̂₁, n̂₂ from PCG64(0) and PCG64(10⁶)).
# As in playground/generate_torus_cluster.py the center of the slice is shifted
# by 6.75 n̂₂ from X₀ (first plane of the torus used for the video).
const R_16 = [1.0062865110546697, 0.9933947568354349, 1.0320211325221642, 1.005245005857652, 0.9732165313419444, 1.0180797527454741, 1.065200002256507, 1.047354048156462, 1.3507385334935105, 1.3114204970267762, 1.3563707876223852, 1.402892818554307, 1.2372478457752816, 1.3846845835247217, 1.3127862336922853, 1.3487412851707583]
const X0_16 = [7.672313765559016, 4.812745383453158, 2.6639710538942083, 3.7570298599183696, 0.25171768811723527, 1.10468440369292, 5.9608046433370765, 5.752504921064107, 5.469813429971451, 3.4102947892722617, 8.863640340845253, 8.718095723367613, 6.093398562236765, 5.781568027239466, 6.119217222582377, 3.4569046088293707, 1.20079713292885, 6.4129055759283125, 4.66958019451339, 2.7575661903604587, 4.318318270808653, 7.906158942336701, 8.302189429733842, 3.1802410159514665, 5.0800084133172705, 2.8609166604082725, 5.2823999573115525, 3.0035035377571577, 3.4808818552091116, 7.913149857060337, 2.0190765630291434, 5.539161331323685]
const N1_16 = [0.027339628754419407, -0.028725774023372848, 0.1392578280446584, 0.022810190217743196, -0.11647956768054132, 0.07862767180943646, 0.2835505801196546, 0.20593968345206415, -0.15302494436938494, -0.2751617943196974, -0.13552901020689365, 0.008986200162864525, -0.5055704004052376, -0.04757553764288904, -0.270919294205836, -0.15922916107577814, -0.11834735044746258, -0.06877844304174131, 0.08950772495732671, 0.22669115065748796, -0.027949445538269075, 0.29713304936888396, -0.14464442410414938, 0.0764347245730572, 0.19645666047704297, 0.020442669206976398, -0.16167150013538584, -0.20042619333115042, -0.09953099609806319, 0.04788071536701551, -0.21953819918249343, -0.04548455026764762]
const N2_16 = [0.08815050094807968, -0.10680864649344908, 0.1319688967427538, 0.11523101283543001, -0.20645098282949867, -0.06939216989933841, 0.4171676987215871, -0.010879947585195593, 0.13395070262773578, -0.30486746961026506, 0.12886889008805977, -0.043747798504101175, 0.0898473481394361, -0.131349775633617, 0.02767449170410305, -0.05535046343203631, -0.17738643240372146, -0.21110872643817546, -0.02085811462950496, -0.10055328767241489, -0.28651747139378725, -0.2984375902661614, -0.011989827626257844, 0.12313747061723948, -0.06095344244651089, -0.30252268186331593, 0.06089286828476715, 0.18970804601709954, 0.1999883187414181, 0.2737841782086356, 0.19115671595352682, -0.10252272073531103]

struct HertzDisks
    R::Vector{Float64}
    L::Float64
end

# Box side L fixed by the packing fraction ϕ.
hertz_disks(R, ϕ) = HertzDisks(R, sqrt(π*sum(abs2, R)/ϕ))

# Minimum image convention.
wrap(d, L) = d - L*round(d/L)

# Steepest descent rule: dX = -∇E.
function (s::HertzDisks)(dX, X, p, t)
    R, L = s.R, s.L
    fill!(dX, 0.)
    N = length(R)
    for i in 1:N, j in i+1:N
        dx = wrap(X[2i-1] - X[2j-1], L)
        dy = wrap(X[2i] - X[2j], L)
        r = sqrt(dx^2 + dy^2)
        σ = R[i] + R[j]
        r < σ || continue
        # V'(r)/r
        f = -2.5/σ*(1 - r/σ)^1.5/r
        dX[2i-1] -= f*dx; dX[2i] -= f*dy
        dX[2j-1] += f*dx; dX[2j] += f*dy
    end
    return nothing
end

# Jacobian of the flow, J = -H with H the Hessian of E.
function hertz_jac!(J, X, s::HertzDisks, t)
    R, L = s.R, s.L
    fill!(J, 0.)
    N = length(R)
    for i in 1:N, j in i+1:N
        dx = wrap(X[2i-1] - X[2j-1], L)
        dy = wrap(X[2i] - X[2j], L)
        r = sqrt(dx^2 + dy^2)
        σ = R[i] + R[j]
        r < σ || continue
        d1 = -2.5/σ*(1 - r/σ)^1.5        # V'(r)
        d2 = 3.75/σ^2*sqrt(1 - r/σ)      # V''(r)
        n = (dx/r, dy/r)
        for a in 1:2, b in 1:2
            k = d2*n[a]*n[b] + d1/r*((a == b) - n[a]*n[b])
            J[2i-2+a, 2i-2+b] -= k; J[2j-2+a, 2j-2+b] -= k
            J[2i-2+a, 2j-2+b] += k; J[2j-2+a, 2i-2+b] += k
        end
    end
    return nothing
end

hertz_energy(X, s::HertzDisks) = sum(
    (r = hypot(wrap(X[2i-1] - X[2j-1], s.L), wrap(X[2i] - X[2j], s.L));
     σ = s.R[i] + s.R[j]; r < σ ? (1 - r/σ)^2.5 : 0.)
    for i in eachindex(s.R) for j in i+1:length(s.R))

# Integrate the flow until the gradient norm drops below gtol. With rtol = 1e-7
# the solver noise can stall ‖∇E‖ slightly above gtol: the final state is then
# accepted if ‖∇E‖ < gacc (still far below the dtol used to match minima).
function relax(X0, s::HertzDisks; rtol = 1e-7, atol = 1e-10, gtol = 1e-10, gacc = 1e-7, T = 1e6)
    g = similar(X0)
    converged(u, t, integ) = (s(g, u, nothing, t); norm(g) < gtol)
    cb = DiscreteCallback(converged, terminate!; save_positions = (false, false))
    prob = ODEProblem(ODEFunction(s; jac = hertz_jac!), copy(X0), (0., T), s)
    sol = solve(prob, CVODE_BDF(); reltol = rtol, abstol = atol, callback = cb,
        save_everystep = false, maxiters = Int(1e7), verbose = false)
    s(g, sol.u[end], nothing, 0.)
    return norm(g) < gacc ? sol.u[end] : nothing
end

# The origin is strictly inside the convex hull of the 2d vectors vs if no
# angular gap between consecutive vectors is larger than π.
function origin_in_hull(vs)
    θ = sort([atan(v[2], v[1]) for v in vs])
    gaps = [diff(θ); 2π + θ[1] - θ[end]]
    return maximum(gaps) < π
end

# Rattlers (as in basinerror/src/julia/check_same_structure.jl): particles with
# less than d+1 = 3 contacts or whose contacts do not surround them, removed
# iteratively. Also returns whether the remaining backbone is jammed.
function find_rattlers(X, s::HertzDisks)
    N = length(s.R)
    contacts = [Int[] for _ in 1:N]
    disp = [NTuple{2,Float64}[] for _ in 1:N]
    for i in 1:N, j in i+1:N
        dx = wrap(X[2j-1] - X[2i-1], s.L)
        dy = wrap(X[2j] - X[2i], s.L)
        if hypot(dx, dy) < s.R[i] + s.R[j]
            push!(contacts[i], j); push!(disp[i], (dx, dy))
            push!(contacts[j], i); push!(disp[j], (-dx, -dy))
        end
    end
    rattler = falses(N)
    changed = true
    while changed
        changed = false
        for i in 1:N
            rattler[i] && continue
            vs = [disp[i][k] for k in eachindex(contacts[i]) if !rattler[contacts[i][k]]]
            if length(vs) < 3 || !origin_in_hull(vs)
                rattler[i] = true; changed = true
            end
        end
    end
    Ns = count(!, rattler)
    nc = sum(count(j -> !rattler[j], contacts[i]) for i in 1:N if !rattler[i]; init = 0)
    jammed = Ns ≥ 2 && nc ≥ 2*(2*(Ns - 1) + 1)
    return rattler, jammed
end

# Minimum stored as its energy, rattler mask and the non-rattler positions
# relative to the first non-rattler particle (removes the translation symmetry).
struct Minimum
    E::Float64
    rattler::BitVector
    pos::Vector{Float64}
end

function Minimum(X, s::HertzDisks, rattler)
    k = findfirst(!, rattler)
    pos = Float64[]
    for i in findall(!, rattler)
        push!(pos, wrap(X[2i-1] - X[2k-1], s.L), wrap(X[2i] - X[2k], s.L))
    end
    return Minimum(hertz_energy(X, s), rattler, pos)
end

function same_minimum(a::Minimum, b::Minimum, L; dtol = 1e-2)
    a.rattler == b.rattler || return false
    for i in 1:2:length(a.pos)
        hypot(wrap(a.pos[i] - b.pos[i], L), wrap(a.pos[i+1] - b.pos[i+1], L)) > dtol && return false
    end
    return true
end

# Thread safe list of the minima found so far. Returns the label of m.
function get_label!(minima, lk, m::Minimum, L)
    lock(lk) do
        for (k, mk) in enumerate(minima)
            abs(mk.E - m.E) < 1e-6*mk.E && same_minimum(mk, m, L) && return k
        end
        push!(minima, m)
        return length(minima)
    end
end

# Label of the minimum reached from X: -1 not converged, 0 fluid state.
function label_minimum(X, s, minima, lk)
    Xf = relax(X, s)
    isnothing(Xf) && return -1
    rattler, jammed = find_rattlers(Xf, s)
    jammed || return 0
    return get_label!(minima, lk, Minimum(Xf, s, rattler), s.L)
end


function compute_soft_spheres(di::Dict)
    @unpack ϕ, res = di
    s = hertz_disks(R_16, ϕ)
    Xc = X0_16 .+ 6.75.*N2_16
    xg = range(-6, 6, length = res)
    yg = range(-3.375, 3.375, length = round(Int, res*6.75/12))
    bsn = zeros(Int, length(xg), length(yg))
    minima = Minimum[]; lk = ReentrantLock()
    p = Progress(length(xg))
    Threads.@threads for i in eachindex(xg)
        for j in eachindex(yg)
            bsn[i,j] = label_minimum(Xc .+ xg[i].*N1_16 .+ yg[j].*N2_16, s, minima, lk)
        end
        next!(p)
    end
    grid = (xg, yg)
    return @strdict(bsn, grid, res, minima)
end


# Paper: 1920 × 1080 slice of 12µs × 6.75µs (µs = 1) for N = 16.
let res = 1920, N = 16, ϕ = 0.9
    params = @strdict res N ϕ
    print_fig(params, "soft_spheres", compute_soft_spheres; cmap = :glasbey_bw_n256,
        xlab = L"p", ylab = L"q", force = false, png_size = (500, 280))
end
