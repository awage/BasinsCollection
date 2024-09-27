using DrWatson
@quickactivate
using Graphs
using LaTeXStrings
using OrdinaryDiffEq:Vern9
using Attractors
using JLD2
include(srcdir("print_fig.jl"))

mutable struct GridParameters{M}
    N::Int
    α::Float64
    incidence::M
    P::Vector{Float64}
    K::Float64
end
# function KuramotoParameters(; N = 10, α = 0.1, K = 6.0, seed = 53867481290)
#     rng = Random.Xoshiro(seed)
#     g = random_regular_graph(N, 3; rng)
#     incidence = incidence_matrix(g, oriented=true)
#     P = [isodd(i) ? +1.0 : -1.0 for i = 1:N]
#     return KuramotoParameters(N, α, incidence, P, K)
# end

function second_order_kuramoto!(du, u, p, t)
    (; N, α, K, incidence, P) = p
    ωs = view(u, N+1:2N)
    du[1:N] .= ωs
    sine_term = K .* (incidence * sin.(incidence' * u[1:N]))
    @. du[N+1:end] .= P - α*ωs - sine_term
    return nothing
end

# function second_order_kuramoto!(du, u, p, t)
#     N = p[1]; α = p[2]; K = p[3]; incidence = p[4]; P = p[5];   
#     du[1:N] .= u[1+N:2*N]
#     du[N+1:end] .= P .- α .* u[1+N:2*N] .- K .* (incidence * sin.(incidence' * u[1:N]))
# end


function get_basins(ni , res)
	# Set up the parameters for the network
	N = 120 # in this case this is the number of oscillators, the system dimension is twice this value
        @load "param_GB_directed.jld2"
        p = GridParameters(N, 0.1, AI, vec(P), K)
        diffeq = (alg = Vern9(), reltol = 1e-6)
        ds = CoupledODEs(second_order_kuramoto!, zeros(2*N), p; diffeq)

        # get the equilibrium state of the system after a long transient.
        uu = trajectory(ds, 1500; Δt = 0.1)
        Δϕ = uu[end][1:N]; Δω = uu[end][N+1:2N]; 
        
        _complete(y) = (length(y) == 5) ? [Δϕ; Δω] : y; 
        # _proj_state(y) = y[N+1:2*N]
        _proj_state(y) = y[N+1:N+5]
        psys = ProjectedDynamicalSystem(ds, _proj_state, _complete)
        # yg = range(-17, 17; length = 31)
        yg = range(-17, 17; length = 2001)
        grid = ntuple(x -> yg, dimension(psys))
	mapper = AttractorsViaRecurrences(psys, grid;  Δt = 1.0, Ttr = 400)
            # consecutive_recurrences = 100,
            # force_non_adaptive = true)
            # consecutive_attractor_steps = 100,
            # Ttr = 600.)

        n = ni + 1; # Indices are 0 based in the paper. 5 is 6 and so on. 
        ϕ = range(-pi, pi; length = res)
        ω = range(-14, 14; length = res)
        ic = [Δϕ; Δω]
        ic_ref = [Δϕ; Δω]
        basins = zeros(Int16,length(ω), length(ϕ))
        for k in 1:length(ω), j in 1:length(ϕ)
           # Prepare the perturbation of the steady state on node n
           # Notice that we are working with a projected system on variables N+1:2N 
           # but we need to define the initial conditions on the full state space. 
           ic[n] = ic_ref[n] + ϕ[j]; ic[n+N] = ic_ref[n+N] + ω[k]
           @show basins[j,k] = mapper(ic)
        end     

	return basins, mapper.bsn_nfo.attractors, (ϕ,ω)
end


function compute_kur_halekotte(di::Dict)
    @unpack ni, res = di
    bsn, att, grid = get_basins(ni, res)
    return @strdict(bsn, grid, ni, res)
end


let res = 300
# ni = 5; 
# params = @strdict ni res
# print_fig(params, "basins_kur", compute_kur_halekotte; xlab = L"\phi", ylab = L"\omega") 
ni = 14
params = @strdict ni res
print_fig(params, "basins_kur", compute_kur_halekotte; xlab = L"\Delta\phi", ylab = L"\omega", force = false) 
# ni = 58
# params = @strdict ni res
# print_fig(params, "basins_kur", compute_kur_halekotte; xlab = L"\phi", ylab = L"\omega") 
# ni = 76
# params = @strdict ni res
# print_fig(params, "basins_kur", compute_kur_halekotte; xlab = L"\phi", ylab = L"\omega") 
end 



