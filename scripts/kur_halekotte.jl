using Graphs
using OrdinaryDiffEq:Tsit5
using DynamicalSystems
using JLD2


function second_order_kuramoto!(du, u, p, t)
    N = p[1]; α = p[2]; K = p[3]; incidence = p[4]; P = p[5];   
    du[1:N] .= u[1+N:2*N]
    du[N+1:end] .= P .- α .* u[1+N:2*N] .- K .* (incidence * sin.(incidence' * u[1:N]))
end


function get_basins(ni , res)

	# Set up the parameters for the network
	N = 120 # in this case this is the number of oscillators, the system dimension is twice this value
        @load "param_GB_directed.jld2"
        ds = ContinuousDynamicalSystem(second_order_kuramoto!, zeros(2*N), [N, 0.1, K, AI, vec(P)], (J,z0, p, n) -> nothing)
        diffeq = (alg = Tsit5(), reltol = 1e-9, maxiters = 1e6)

        # get the equilibrium state of the system after a long transient.
        uu = trajectory(ds, 1500; Δt = 0.1, diffeq)
        Δϕ = uu[end][1:N]; Δω = uu[end][N+1:2N]; 
        
        _complete(y) = (length(y) == N) ? [Δϕ; Δω] : y; 
        _proj_state(y) = y[N+1:2*N]
        psys = projected_integrator(ds, _proj_state, _complete; diffeq)
        yg = range(-17, 17; length = 31)
        grid = ntuple(x -> yg, dimension(psys))
	mapper = AttractorsViaRecurrences(psys, grid; sparse = true, Δt = .1,   
            mx_chk_fnd_att = 100,
            mx_chk_loc_att = 100,
            Ttr = 600.)

        n = ni + 1; # Indices are 0 based in the paper. 5 is 6 and so on. 
        ϕ = range(-pi, pi; length = res)
        ω = range(-14, 14; length = res)
        ic = [Δϕ; Δω]
        basins = zeros(Int16,length(ω), length(ϕ))
        for k in 1:length(ω), j in 1:length(ϕ)
           # Prepare the perturbation of the steady state on node n
           # Notice that we are working with a projected system on variables N+1:2N 
           # but we need to define the initial conditions on the full state space. 
           ic[n] = ϕ[j]; ic[n+N] = ω[k]
           @show basins[j,k] = mapper(ic)
        end     

	return basins, mapper.bsn_nfo.attractors, (ϕ,ω)
end


function get_trajectory(ui)
	# Set up the parameters for the network
	N = 120 # in this case this is the number of oscillators, the system dimension is twice this value
        @load "param_GB_directed.jld2"
        ds = ContinuousDynamicalSystem(second_order_kuramoto!, zeros(2*N), [N, α, K, AI, vec(P)], (J,z0, p, n) -> nothing)
        diffeq = (alg = Vern9(), abstol = 1e-12, reltol = 1e-12,  maxiters = 1e7)
        pp = trajectory(ds, 1500, ui; Δt = 0.1, diffeq)
	return pp
end

