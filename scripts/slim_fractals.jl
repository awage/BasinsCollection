using DrWatson
@quickactivate
using DynamicalSystems
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter

function slim_fractal_U3!(dz, z, p, n)
    x = z[1]; dx = z[2]; 
    y = z[3]; dy = z[4]; 
    μ = p[1] 
    # U(r,θ) = -r^2*cos(3θ) + 1/2*r^4
    dUdx = 2*x^3+(-2*x^4-3*x^2*y^2+3*y^4)/(x^2+y^2)^(3/2) + 2*x*y^2
    dUdy = 2*y^3 + y*(x*(7*x^2+3*y^2)/(x^2+y^2)^(3/2) + 2*x^2)
    # dUdx = 2*x
    # dUdy = 2*y
    dz[1] = dx
    dz[2] = -μ*dx - dUdx
    dz[3] = dy 
    dz[4] = -μ*dy - dUdy
    return 
end



function compute_slim_fractal(di)
    @unpack res,μ = di
    p = [μ]
    df = ContinuousDynamicalSystem(slim_fractal_U3!, rand(4), p) 
    x1 = x2 = y1 = y2 =  range(-2, 2, length = 1001)
    grid_rec = (x1, x2, y1, y2)
    diffeq = (reltol = 1e-9,  alg = Vern9())
    mapper = AttractorsViaRecurrences(df, grid_rec; Δt = 1.,
            mx_chk_fnd_att = 300, 
            mx_chk_loc_att = 300, 
            sparse = true, diffeq
    )
    x1 =  y1 = range(-2, 2, length = res)
    grid = (x1, y1)

    bsn = @showprogress [ mapper([x,0,y,0.]) for x in x1, y in y1]
    att = mapper.bsn_nfo.attractors

    return @strdict(bsn, att, grid, μ, res)
end



function print_fig(w, h, cmap, μ, res) 
    params = @strdict μ res

    data, file = produce_or_load(
        datadir("basins"), params, compute_slim_fractal;
        prefix = "slim_fractal", storepatch = false, suffix = "jld2", force = false
    )


    @unpack bsn, grid = data
    xg ,yg = grid
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"\dot{\theta}", xlabel = L"\theta", yticklabelsize = 30, 
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
    save(string("../plots/slim_fractal_", res, ".png"),fig)
end

μ = 0.2
print_fig(600,600, nothing, μ, 100) 
# using Plots
# df = ContinuousDynamicalSystem(slim_fractal_U3!, rand(4), [μ]) 
# diffeq = (reltol = 1e-9,  alg = Vern9())
# u = trajectory(df, 100, [0.4, 0, 0.4, 0]; Ttr = 100, diffeq)
# plot(Matrix(u[1000:2000])[:,4])
