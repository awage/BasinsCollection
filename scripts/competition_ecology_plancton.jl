using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter

function comptetion_model!(du, u, p, t)
    n = length(u)-3
    R = u[n+1:n+3]
    d = 1.
    D = 0.25/d
    r = 1/d
    K = [0.2 0.05 0.50 0.05 0.50 0.03 0.51 0.51;
        0.15 0.06 0.05 0.50 0.30 0.18 0.04 0.31 ;
        0.15 0.50 0.30 0.06 0.05 0.18 0.31 0.04]
    C = [0.2  0.10 0.10 0.10 0.10 0.22 0.10 0.10;
        0.10 0.20 0.10 0.10 0.20 0.10 0.22 0.10;
        0.10 0.10 0.20 0.20 0.10 0.10 0.10 0.22]
    tmp_v = zeros(n)
    for k in 1:n
        tmp_v[k] = r*u[k]*min(R[1]/(K[1,k]+R[1]), R[2]/(K[2,k]+R[2]), R[3]/(K[3,k]+R[3]))
        du[k] = tmp_v[k] - u[k]*D
    end
    du[n+1:n+3] = D*(10 .- R) .- C*tmp_v
end

# Reference: Fundamental unpredictability in multispecies competition
# Huisman, J; Weissing, FJ, American Naturalist DOI: 10.1086/319929
function compute_competition_model(di)
    @unpack res = di
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    ds = CoupledODEs(comptetion_model!, rand(8+3), []; diffeq)
    yg = range(0, 100; length = 41)
    grid = ntuple(x -> yg, dimension(ds))
    mapper = AttractorsViaRecurrences(ds, grid;
            mx_chk_fnd_att = 300, 
            mx_chk_loc_att = 300, 
    )
    xg = range(0.,2.,length = res)
    yg = range(0.,2.,length = res)
    bsn = zeros(Int16, res, res)
# Caching initial conditions with parallel threads. This is the longest part.
    ds_arr = [deepcopy(ds) for k in 1:Threads.nthreads()]
    ics = Array{Vector{Float64}}(undef, res, res)
    @showprogress for i in 1:res
        Threads.@threads for j in 1:res
            u,t = trajectory(ds_arr[Threads.threadid()], 1000., [0.1, xg[i], 0.1, yg[j], 0.1, 0., 0., 0., 10., 10., 10.])
            ics[i,j] = vec(u[end]) + [zeros(5); 0.1; 0.1; 0.1; zeros(3)]
        end
    end
    bsn = @showprogress  [mapper(ics[i,j]) for i in 1:res, j in 1:res]
    att = mapper.bsn_nfo.attractors
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end

let res = 1200
params = @strdict res
print_fig(params, "comptetion_model", compute_competition_model; xlab = L"N_2", ylab = L"N_4", force = false)
end
