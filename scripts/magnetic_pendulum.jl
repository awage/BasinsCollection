using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq

# https://cdn.ima.org.uk/wp/wp-content/uploads/2020/03/Chaos-in-the-Magnetic-Pendulum-from-MT-April-2020.pdf
#ds = mag_pendulum(γ=1, d=0.5, α=0.175, ω=1., N=4)
struct MagneticPendulum
    magnets::Vector{SVector{2, Float64}}
end
mutable struct MagneticPendulumParams
    γs::Vector{Float64}
    d::Float64
    α::Float64
    ω::Float64
end

function (m::MagneticPendulum)(u, p, t)
    x, y, vx, vy = u
    γs::Vector{Float64}, d::Float64, α::Float64, ω::Float64 = p.γs, p.d, p.α, p.ω
    dx, dy = vx, vy
    dvx, dvy = @. -ω^2*(x, y) - α*(vx, vy)
    for (i, ma) in enumerate(m.magnets)
        δx, δy = (x - ma[1]), (y - ma[2])
        D = sqrt(δx^2 + δy^2 + d^2)
        dvx -= γs[i]*(x - ma[1])/D^3
        dvy -= γs[i]*(y - ma[2])/D^3
    end
    return SVector(dx, dy, dvx, dvy)
end

"""
    magnetic_pendulum(u=[0.7,0.7,0,0]; d=0.3, α=0.2, ω=0.5, N=3, γs=fill(1.0,N))

Create a pangetic pendulum with `N` magnetics, equally distributed along the unit circle,
with dynamical rule
```math
\\begin{aligned}
\\ddot{x} &= -\\omega ^2x - \\alpha \\dot{x} - \\sum_{i=1}^N \\frac{\\gamma_i (x - x_i)}{D_i^3} \\\\
\\ddot{y} &= -\\omega ^2y - \\alpha \\dot{y} - \\sum_{i=1}^N \\frac{\\gamma_i (y - y_i)}{D_i^3} \\\\
D_i &= \\sqrt{(x-x_i)^2  + (y-y_i)^2 + d^2}
\\end{aligned}
```
where α is friction, ω is eigenfrequency, d is distance of pendulum from the magnet's plane
and γ is the magnetic strength.
"""
function magnetic_pendulum(u = [sincos(0.12553*2π)..., 0, 0];
    γ = 1.0, d = 0.3, α = 0.2, ω = 0.5, N = 3, γs = fill(γ, N), diffeq)
    m = MagneticPendulum([SVector(cos(2π*i/N), sin(2π*i/N)) for i in 1:N])
    p = MagneticPendulumParams(γs, d, α, ω)
    return CoupledODEs(m, u, p; diffeq)
end


function compute_mag_pend(di::Dict)
    @unpack γ, d, α, ω, N, res = di
    diffeq = (alg = Vern9(), reltol = 1e-6, maxiters = 1e8)
    ds = magnetic_pendulum(γ=γ, d=d, α=α, ω=ω, N=N, diffeq = diffeq)
    xg = yg = range(-10.,10.,length = 5001)
    psys = ProjectedDynamicalSystem(ds, [1,2], [0., 0.])
    mapper = AttractorsViaRecurrences(psys, (xg, yg) ; Δt = 0.1)
    xg = yg = range(-5.,5.,length = res)
    bsn, att = basins_of_attraction(mapper, (xg,yg))
    grid = (xg,yg)
    return @strdict(bsn, att, grid, γ, d, α, ω, N, res)
end



γ=1; d=0.3; α=0.2; ω=0.5; N=3; #res = 300
params = @strdict γ d α ω N res
print_fig(params, "mag_pend", compute_mag_pend; ylab = L"y", xlab = L"x", force = false)
