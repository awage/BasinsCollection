using DrWatson
@quickactivate "PerspectiveFigures" # exports DynamicalSystems, GLMakie and other goodies in `src`
using DynamicalSystems
using CairoMakie
using LaTeXStrings

function newton_map(dz,z, p, n)
    f(x) = x^p[1]-1
    df(x)= p[1]*x^(p[1]-1)
    z1 = z[1] + im*z[2]
    dz1 = f(z1)/df(z1)
    z1 = z1 - dz1
    dz[1]=real(z1)
    dz[2]=imag(z1)
    return
end

# dummy function to keep the initializator happy
function newton_map_J(J,z0, p, n)
   return
end

function compute_basins_newton(di::Dict)
    @unpack N, res = di
    ds = DiscreteDynamicalSystem(newton_map,[0.1, 0.2], [N] , newton_map_J)
    xg=range(-3.,3.,length = res)
    yg=range(-3.,3.,length = res)
    mapper = AttractorsViaRecurrences(ds, (xg, yg))
    bsn, att = basins_of_attraction(mapper)
    return @strdict(bsn, att, (xg,yg), N, res)
end


function print_fig(w,h,cmap, N, res)
    params = @strdict N res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_newton;
        prefix = "newton", storepatch = false, suffix = "jld2", force = false
    )

    @unpack bsn, grid = data
# Remove spurious point 
    ind = findall(bsn .== -1)
    bsn[ind] .= 1
    xg, yg = grid
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"\Im{(z)}", xlabel = L"\Re{(z)}", 
            aspect = AxisAspect(1.),
            yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
    save("../plots/basins_newton.svg",fig)

end


function get_Sb(N, res)
    params = @strdict N res
    data, file = produce_or_load(
        datadir("basins"), params, compute_basins_newton;
        prefix = "newton", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    ind = findall(bsn .== -1)
    bsn[ind] .= 1
    @show  basin_entropy(bsn)
end

cmap = ColorScheme([RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )
print_fig(600, 600, cmap, 3, 2000) 
