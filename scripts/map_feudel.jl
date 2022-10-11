using DrWatson
@quickactivate
using DynamicalSystems
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9



# Basin bifurcation in quasiperiodically forced systems Ulrike Feudel, Annette Witt, Ying-Cheng Lai, and Celso Grebogi PRE 28, 1998
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

# dummy function to keep the initializator happy
function chaotic_map_J(J, z0, p, n)
    return
end



function compute_feudel(di::Dict)
    @unpack a, r, res = di
    ds = DiscreteDynamicalSystem(chaotic_map, [1.0, 0.0], [a, r], chaotic_map_J)
    θ = range(0.0, 1.0, length = 2500)
    xg = range(0.0, 1.0, length = 2500)
    mapper = AttractorsViaRecurrences(ds, (xg,θ); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, safety_counter_max = Int(1e7), show_progress = true)
    θ = range(0.0, 1.0, length = res)
    xg = range(0.0, 1.0, length = res)
    bsn, att = basins_of_attraction(mapper, (xg,θ); show_progress = true)
    grid = (xg, θ)
    return @strdict(bsn, att, grid, b, res)
end


function print_fig(w, h, cmap, a, r, res)
    params = @strdict res a r
    data, file = produce_or_load(
        datadir("basins"), params, compute_feudel;
        prefix = "feudel", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    xg, yg = grid
    fig = Figure(resolution = (w, h))
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
    save(string(projectdir(), "/plots/feudel",res,".png"),fig)
end

r = 3.833
a = 0.0015
print_fig(600,600, nothing, a, r, 400) 
