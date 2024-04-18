using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie


# https://cdn.ima.org.uk/wp/wp-content/uploads/2020/03/Chaos-in-the-Magnetic-Pendulum-from-MT-April-2020.pdf
#ds = mag_pendulum(γ=1, d=0.5, α=0.175, ω=1., N=4)

function compute_mag_pend(di::Dict)
    @unpack γ, d, α, ω, N, res = di
    ds = Systems.magnetic_pendulum(γ=γ, d=d, α=α, ω=ω, N=N)
    xg = yg = range(-2.,2.,length = res)
    psys = ProjectedDynamicalSystem(ds, [1,2], [0., 0.])
    mapper = AttractorsViaRecurrences(psys, (xg, yg); Δt = 0.1,  mx_chk_fnd_att = 1000, mx_chk_att = 10, mx_chk_hit_bas = 100)
    bsn, att = basins_of_attraction(mapper)
    grid = (xg,yg)
    return @strdict(bsn, att, grid, γ, d, α, ω, N, res)
end



γ=1; d=0.3; α=0.2; ω=0.5; N=3; res = 300
params = @strdict γ d α ω N res
print_fig(params, "mag_pend", compute_mag_pend; ylab = L"\dot{x}", xlab = L"x")
