using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))

function slim_fractal_U2!(dz, z, p, n)
    x = z[1]; dx = z[2]; 
    y = z[3]; dy = z[4]; 
    μ = p[1] 
    # U(r,θ) = -r^3*cos(3θ) + 3/4*r^4
    # U(x,y) = -xr^2(-3 + 4x^2/r^2) + 3/4*r^4
    dUdx(x,y) = 3*(x^3-x^2+x*y^2+y^2)
    dUdy(x,y) = 3*y*(x^2+2*x+y^2)
    dz[1] = dx
    dz[2] = -μ*dx - dUdx(x,y)
    dz[3] = dy 
    dz[4] = -μ*dy - dUdy(x,y)
    return 
end


function compute_slim_fractal(di)
    @unpack res,μ = di
    p = [μ]
    diffeq = (reltol = 1e-9,  alg = Vern9())
    df = CoupledODEs(slim_fractal_U2!, rand(4), p; diffeq) 
    x1 = x2 = y1 = y2 =  range(-2, 2, length = 1001)
    grid_rec = (x1, x2, y1, y2)
    mapper = AttractorsViaRecurrences(df, grid_rec; Δt = 1.,
            mx_chk_fnd_att = 300, 
            mx_chk_loc_att = 300, 
    )
    θr = range(0, 2π/3, length = res)
    ur = range(0, 1, length = res)
    grid = (θr, ur)
    # vmax for potential U₂
    # ∂U∂r = -3r²*cos(3θ) + 3*r^3
    vmax(θ) = sqrt(2*(-3*2^2*cos(3*θ) + 3*2^3))
    u0(r,θ) = [2cos(θ), vmax(θ)*r*cos(θ+π/2), 2sin(θ),  vmax(θ)*r*sin(θ+π/2)]
    bsn = @showprogress [ mapper(u0(u,θ))  for θ in θr, u in ur ]
    att = mapper.bsn_nfo.attractors
    return @strdict(bsn, att, grid, μ, res)
end

μ = 0.2; #res = 800
params = @strdict μ res
print_fig(params, "slim_fractal", compute_slim_fractal; ylab = L"v_0/v_{max}", xlab = L"\theta_0/(2\pi/3)")
