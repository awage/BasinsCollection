using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter

# function comptetion_model!(du, u, p, t)

    # K = [0.2 0.05 1.00 0.05 1.20;
    #      0.25 0.10 0.05 1.00 0.40;
    #      0.15 0.95 0.35 0.10 0.05]

    # C = [0.2 0.10 0.10 0.10 0.10;
    #     0.10 0.20 0.10 0.10 0.20;
    #     0.10 0.10 0.20 0.20 0.1]
    K = [0.2 0.05 0.50 0.05 0.50 0.03 0.51 0.51;
        0.15 0.06 0.05 0.50 0.30 0.18 0.04 0.31 ;
        0.15 0.50 0.30 0.06 0.05 0.18 0.31 0.04]

    C = [0.2  0.10 0.10 0.10 0.10 0.22 0.10 0.10;
        0.10 0.20 0.10 0.10 0.20 0.10 0.22 0.10;
        0.10 0.10 0.20 0.20 0.10 0.10 0.10 0.22]
    d = 1.
    D = 0.25/d
    S = 10. *ones(3)
    n = 8 
    r = ones(n)*d
    m = D*ones(n)
# array of functions
    μ = [ R -> min(r[k]*R[1]/(K[1,k]+R[1]), r[k]*R[2]/(K[2,k]+R[2]), r[k]*R[3]/(K[3,k]+R[3]))   for k in 1:n]


function comptetion_model!(du, u, p, t)
    n = length(u)-3
    R = u[n+1:n+3]
    # μ_v = [min(r[k]*R[1]/(K[1,k]+R[1]), r[k]*R[2]/(K[2,k]+R[2]), r[k]*R[3]/(K[3,k]+R[3]))   for k in 1:n]
    μ_v = [min(R[1]/(K[1,k]+R[1]), R[2]/(K[2,k]+R[2]), R[3]/(K[3,k]+R[3]))   for k in 1:n]

    tmp_v = u[1:n].*μ_v[1:n]
    du[1:n] = tmp_v .- u[1:n]*D
    # for k in 1:n
    #     du[k] = u[k]*(μ_v[k] - m[k])
    # end

    tv = C*tmp_v    
    du[n+1:n+3] = D*(10 .- R) .- tv
    # for k in 1:3
    #     du[n+k] = D*(10 - R[k]) - tv[k] # sum(u[1:n].*μ_v[1:n].*C[k,1:n])
    # end

end

# Reference: Fundamental unpredictability in multispecies competition
# Huisman, J; Weissing, FJ, American Naturalist DOI: 10.1086/319929
function compute_competition_model(di)
    @unpack res = di
    ds = ContinuousDynamicalSystem(comptetion_model!, rand(8+3), [])
    yg = range(0, 100; length = 41)
    grid = ntuple(x -> yg, dimension(ds))
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    mapper = AttractorsViaRecurrences(ds, grid;
            # mx_chk_lost = 10, 
            mx_chk_fnd_att = 300, 
            mx_chk_loc_att = 300, 
            # mx_chk_att = 2,
             sparse = true, diffeq
    )
    xg = range(0.,2.,length = res)
    yg = range(0.,2.,length = res)
    bsn = zeros(Int16, res, res)
# Caching initial conditions with parallel threads. This is the longest part.
    ics = Array{Vector{Float64}}(undef, res, res)
    for i in 1:res
        Threads.@threads for j in 1:res
            u = trajectory(ds, 1000, [0.1, xg[i], 0.1, yg[j], 0.1, 0., 0., 0., 10., 10., 10.]; diffeq)
            @show ics[i,j] = vec(u[end]) + [zeros(5); 0.1; 0.1; 0.1; zeros(3)]
        end
    end
    # u0 = vec(u[end]) + [zeros(5); 0.1; 0.1; 0.1; zeros(3)]
    bsn = @showprogress  [mapper(ics[i,j]) for i in 1:res, j in 1:res]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end


function print_fig(w, h, cmap, res)
    params = @strdict res
    data, file = produce_or_load(
        datadir("basins"), params, compute_competition_model;
        prefix = "comptetion_model", storepatch = false, suffix = "jld2", force = false
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
    save(string(projectdir(), "/plots/comptetion_model",res,".png"),fig)
end

print_fig(600,600, nothing, 800) 
# diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
# ds = ContinuousDynamicalSystem(comptetion_model!, rand(5+3), [])
# integ = integrator(ds,  [0.1, 0.1, 0.1, 0.1, 0.1, 10., 10., 10.]; diffeq)
# @profview u = trajectory(ds, 10, [0.1, 0.1, 0.1, 0.1, 0.1, 10., 10., 10.]; diffeq)
# @profview step!(integ, 30)
# using Plots
# plot(Matrix(u)[:,1:5])
