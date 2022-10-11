using DrWatson
@quickactivate
using DynamicalSystems
using LaTeXStrings
using CairoMakie

# C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983
function grebogi_map(dz,z, p, n)
    θ = z[1]; x = z[2]
    J₀=0.3; a=1.32; b=0.9;
    dz[1]= θ + a*sin(2*θ) - b*sin(4*θ) -x*sin(θ)
    dz[1] = mod(dz[1],2π) # to avoid problems with attracto at θ=π
    dz[2]=-J₀*cos(θ)
    return
end

# dummy function to keep the initializator happy
function grebogi_map_J(J,z0, p, n)
   return
end





function compute_grebogi(di::Dict)
    @unpack res = di
    ds = DiscreteDynamicalSystem(grebogi_map,[1., -1.], [] , grebogi_map_J)
    θ = range(0, 2π, length = 2000)
    xg = range(-0.5, 0.5, length = 2000)
    mapper = AttractorsViaRecurrences(ds, (θ,xg); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, safety_counter_max = Int(1e7), show_progress = true)
    θ = range(0, 2π, length = res)
    xg = range(-0.5, 0.5, length = res)
    bsn, att = basins_of_attraction(mapper, (θ,xg); show_progress = true)
    grid = (θ,xg)
    return @strdict(bsn, att, grid, b, res)
end


function print_fig(w, h, cmap, res)
    params = @strdict res 
    data, file = produce_or_load(
        datadir("basins"), params, compute_grebogi;
        prefix = "grebogi", storepatch = false, suffix = "jld2", force = false
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
    save(string(projectdir(), "/plots/grebogi",res,".png"),fig)
end

print_fig(600,600, nothing, 400) 
