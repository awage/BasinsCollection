using DrWatson
@quickactivate
using CairoMakie
using LaTeXStrings
using DynamicalSystems
using OrdinaryDiffEq

# Multi dimenioned intertwined basin boundaries: basin structure of
# THE KICKED DOUBLE ROTOR
# Celso GREBOGI a, Eric KOSTELICH b, Edward OTT a,c and James A. YORKE b
# Physica 25D 1987
function kicked_ott!(dz, z, p, t)
    f0 = p[1]; T = 1
    m1 = m2 = 1; 
    ν1 = ν2 = ν = 0.1
    l1 = 1/√2; l2 =1 
    s1 = ν*(3+√5)/2; s2 = ν*(3-√5)/2; 

    u1 = [(ν - s1) , ν]/sqrt(ν^2 + (ν - s1)^2)
    u2 = [(ν - s2) , ν]/sqrt(ν^2 + (ν - s2)^2)

    W1 = u1*u1'
    W2 = u2*u2'

    L = W1*exp(-s1) + W2*exp(-s2)
    M = W1*exp(-s1)/s1 + W2*exp(-s2)/s2
    I1 = (m1 + m2)*l1^2
    I2 = m2*l2^2
    I1 = (m1 + m2)*l1^2; I2 = m2*l2^2

    dz[1:2] = rem2pi.(M*z[3:4] .+ z[1:2], RoundDown)
    dz[3:4] = L*z[3:4] .+ [l1*f0/I1; l2*f0/I2].*sin.(dz[1:2])

    return 
end


function compute_kicked_rotor_ott(di)
    @unpack f0, res = di
    u0 = [0.; 0.; 0.; 0.]
    p = [f0]
    df = DiscreteDynamicalSystem(kicked_ott!, u0, p) 
    θ1 = θ2= range(-2π, 2π, length = 10001)
    θd1 = θd2 = range(-30, 30, length = 10001)
    grid_rec = (θ1, θ2, θd1, θd2)
    mapper = AttractorsViaRecurrences(df, grid_rec;
            # mx_chk_lost = 10, 
            mx_chk_fnd_att = 3000, 
            mx_chk_loc_att = 3000, 
            # mx_chk_att = 2,
             sparse = true,
    )
    θ1 = range(-π, π, length = res)
    θ2 = range(-π, π, length = res)
    bsn = [ mapper([x,y,0., 0.]) for x in θ1, y in θ2]
    att = mapper.bsn_nfo.attractors
    grid = (θ1, θ2)
    return @strdict(bsn, att, grid, f0, res)
end



function print_fig(w, h, cmap, f0, res) 
    params = @strdict f0 res

    data, file = produce_or_load(
        datadir("basins"), params, compute_kicked_rotor_ott;
        prefix = "kicked_rotor_ott", storepatch = false, suffix = "jld2", force = true
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
    save(string("../plots/kicked_rotor_cont_", res, ".png"),fig)
end

 
f0 = 0.1
print_fig(600,600, nothing, f0, 700) 
