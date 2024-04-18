using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
# https://cdn.ima.org.uk/wp/wp-content/uploads/2020/03/Chaos-in-the-Magnetic-Pendulum-from-MT-April-2020.pdf
#ds = mag_pendulum(γ=1, d=0.5, α=0.175, ω=1., N=4)

function compute_mag_pend(di::Dict)
    @unpack γ, d, α, ω, N, res = di
    ds = Systems.magnetic_pendulum(γ=γ, d=d, α=α, ω=ω, N=N)
    xg = yg = range(-2.,2.,length = res)
    psys = projected_integrator(ds, [1,2], [0., 0.])
    diffeq = (;reltol = 1e-9, alg = Vern9(), maxiters = 1e6)
    mapper = AttractorsViaRecurrences(psys, (xg, yg); Δt = 0.1, diffeq, mx_chk_fnd_att = 1000, mx_chk_att = 10, mx_chk_hit_bas = 100)
    bsn, att = basins_of_attraction(mapper)
    grid = (xg,yg)
    return @strdict(bsn, att, grid, γ, d, α, ω, N, res)
end


function print_fig(w,h,cmap, γ, d, α, ω, N, res)
    params = @strdict γ d α ω N res
    data, file = produce_or_load(
        datadir("basins"), params, compute_mag_pend;
        prefix = "mag_pend", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    xg, yg = grid

    fig = Figure(size = (w, h))
    ax = Axis(fig[1,1], ylabel = L"$\dot{x}$", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    if isnothing(cmap)
        heatmap!(ax, xg, yg, bsn, rasterize = 1)
    else
        heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    end
    save(string("../plots/mag_pend", res,".png"),fig)
end

γ=1; d=0.3; α=0.2; ω=0.5; N=3;
print_fig(600, 600, nothing, γ, d, α, ω, N, 300)
