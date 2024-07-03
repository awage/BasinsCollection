using DrWatson
@quickactivate 
using CairoMakie
using LaTeXStrings
using Attractors
using OrdinaryDiffEq:Vern9

# Reference: PHYSICAL REVIEW E 85, 035202(R) (2012)
#How to obtain extreme multistability in coupled dynamical systems
#C. R. Hens,1 R. Banerjee,1,2 U. Feudel,3 and S. K. Dana1
function coupled_rossler!(du, u, p, t)
    a = p[1]; b = p[2]; c = p[3];
    du[1] = - u[2] - u[3]
    du[2] = u[1] + a*u[5]
    du[3] = b - c*u[3] + u[4]*u[6]
    du[4] = u[1] - u[4] -u[2] - u[3]
    du[5] = u[4] + a*u[5]
    du[6] = b + u[6]*(u[4] -c)
end


function compute_EM_rossler(di::Dict)
    @unpack res = di
    a = 0.2; b = 0.2; c = 5.7
    ds = CoupledODEs(coupled_rossler!, zeros(6), [a, b, c], (J,z0, p, n) -> nothing)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    yg = range(-25, 25; length = 10001)
    grid = ntuple(x -> yg, dimension(ds))
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, Δt = .1,   
        mx_chk_fnd_att = 300,
        mx_chk_loc_att = 100, maximum_iterations = Int(1e7), diffeq)

    y1 = range(-8, 8, length = res)
    y2 = 1.
    ics = [ [-0.1, 0.01, 0.3, t1, y2, 2.] for t1 in y1]
    bas = [ mapper(u) for u in ics]
    att = mapper.bsn_nfo.attractors
    
    # Plot "bifurcation" diagram as a function of the variable y1: 
    pnt_lst = Vector{Vector{Float64}}(undef,1)
    for k in 1:length(y1)
        p = mapper.bsn_nfo.attractors[k]
        tra = trajectory(ds, 20, p[1]; Ttr = 20, Δt = 0.1, diffeq)
        for y in tra
             v = [y1[k], y[1]]
             push!(pnt_lst, v)
        end
    end
    P = Dataset(pnt_lst[2:end])
    return @strdict(bas, att, P, y1, res)
end


function print_fig(w,h,res)
    params = @strdict res
    data, file = produce_or_load(
        datadir("bif_diag"), params, compute_EM_rossler;
        prefix = "rossler", storepatch = false, suffix = "jld2", force = false
    )
    @unpack P, bas, y1 = data
    fig = Figure(size = (w, h))
    ax = Axis(fig[1,1], ylabel = "xn", yticklabelsize = 20, xticklabelsize = 20, ylabelsize = 20)
    scatter!(ax, P[:,1],P[:,2], markersize = 0.7, color = :black, rasterize = 4)
    save(string("../plots/bif_diag_em_rossler", res,".png"),fig)

end


print_fig(800,400, 20)
