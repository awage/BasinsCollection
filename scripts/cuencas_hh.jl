using DrWatson
@quickactivate "PerspectiveFigures" # exports DynamicalSystems, GLMakie and other goodies in `src`

using OrdinaryDiffEq
using DynamicalSystems
using CairoMakie
using LaTeXStrings

function henon_heiles_de(dv,v,p,t)
    x, px, y, py = v
    ∂H∂x = x+2*x*y;
    ∂H∂y = y-y^2+x^2;
    ∂H∂px = px;
    ∂H∂py = py;
    dv[1] = ∂H∂px;
    dv[2] = -∂H∂x;
    dv[3] = ∂H∂py;
    dv[4] = -∂H∂y;
end


function salida(sol)
    a=sol[1,end]
    b=sol[3,end]
    if b>1
        sal=1;
    elseif b<0 && a<-1
        sal=2;
    elseif b<0 && a>1
        sal=3;
    else sal=0;
    end
    return sal
end


function get_hh_exit(Ei, y, x)
    Esqrt = 2*Ei -x^2 -2*x^2*y - y^2 + 2/3 * y^3 
    if Esqrt < 0
        return 0
    end
    px = -sqrt(Esqrt)*y/sqrt(x^2+y^2)
    py = sqrt(Esqrt)*x/sqrt(x^2+y^2)
    vi = [x, px, y, py];
    tspan=(0.0,4000)
    # Define callback to halt solver
    condition(u,t,integrator) = (u[1]^2+u[3]^2)>100
    affect!(integrator) = terminate!(integrator)
    cb = DiscreteCallback(condition,affect!)
    prob = ODEProblem(henon_heiles_de, vi, tspan, callback=cb)
    sol = solve(prob, Vern9(), reltol=1e-8, abstol=1e-8)
    return salida(sol)
end


function comp_basins_hh(Ei, res)
    ry = range(-1.1, 1.4, length = res)
    rx = range(-1.2, 1.2, length = res)
    @time basins = [get_hh_exit(Ei, y, x) for x in rx, y in ry]
    return basins, (rx, ry)
end


function makesim(di::Dict)
    @unpack E, res = di
    bsn, grid = comp_basins_hh(E, res)
    return @strdict(bsn, grid, E, res)
end


function print_fig(w, h, E, res) 
    params = @strdict E res
    data, file = produce_or_load(
        datadir("basins"), params, makesim;
        prefix = "hh", storepatch = false, suffix = "jld2", force = false
    )
    @unpack bsn, grid = data
    xg, yg = grid
    cmap = ColorScheme([RGB(1,1,1), RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"y", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    heatmap!(ax, xg, yg, bsn; rasterize = 1, colormap = cmap)
    save(string("../plots/basins_hh_", E, ".png"),fig)

end
