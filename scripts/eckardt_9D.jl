using DrWatson
@quickactivate
using LaTeXStrings
# using Attractors
using DynamicalSystems
using OrdinaryDiffEq:Vern9
using CairoMakie
using Random


mutable struct E9DParameters{M}
    k::M
    σ::M
    Re::Float64
end
function E9DParameters(; Re = 307.)
   Lx = 1.75π; Lz = 1.2π
   α = 2π/Lx; β = π/2; γ = 2π/Lz; 
    Kαγ = sqrt(α^2 + γ^2); 
    Kβγ = sqrt(β^2 + γ^2); 
    Kαβγ = sqrt(α^2 + β^2 + γ^2)
    k = [   β^2; 
            4*β^2/3+ γ^2;
            β^2+γ^2; 
            (3*α^2+4*β^2)/3; 
            α^2 + β^2; 
            (3*α^2 + 4*β^2 + 3*γ^2)/3;  
            α^2 + β^2 + γ^2; 
            α^2 + β^2 + γ^2; 
            9*β^2]
    σ = [-√(3/2)*β*γ/Kαβγ;  √(3/2)*β*γ/Kβγ;
         5√2*γ^2/(3√3*Kαγ); -γ^2/(√6*Kαγ); -α*β*γ/(√6*Kαγ*Kαβγ); -√(3/2)*β*γ/Kβγ; -√(3/2)*β*γ/Kβγ; 
         2*α*β*γ/(√6*Kαγ*Kβγ); (β^2*(3*α^2+γ^2)-3*γ^2*(α^2+γ^2))/(√6*Kαγ*Kβγ*Kαβγ);
         -α/√6; -10*α^2/(3*√6*Kαγ); -√(3/2)*α*β*γ/(Kαγ*Kβγ); -√(3/2)*α^2*β^2/(Kαγ*Kβγ*Kαβγ); -α/√6; 
         α/√6; α^2/(√6*Kαγ); -α*β*γ/(√6*Kαγ*Kαβγ); α/√6; 2*α*β*γ/(√6*Kαγ*Kβγ);
         α/√6; √(3/2)*β*γ/Kαβγ; 10*(α^2 - γ^2)/(3√6*Kαγ); -2√2*α*β*γ/(√3*Kαγ*Kβγ); α/√6; √(3/2)*β*γ/Kαβγ; 
         -α/√6; (γ^2-α^2)/(√6*Kαγ); α*β*γ/(√6*Kαγ*Kβγ);
         2*α*β*γ/(√6*Kαγ*Kαβγ); γ^2*(3*α^2-β^2+3*γ^2)/(√6*Kαγ*Kβγ*Kαβγ);
        √(3/2)*β*γ/Kβγ;  -√(3/2)*β*γ/Kαβγ 
        ] 
    return E9DParameters(k, σ, Re)
end

function E9D!(du, u, p, t)
    (; k, σ, Re) = p
    du[1] = -u[1]*k[1]/Re + σ[1]*u[6]*u[8] + σ[2]*u[2]*u[3] + k[1]/Re; 
    du[2] = -u[2]*k[2]/Re + σ[3]*u[4]*u[6] + σ[4]*u[5]*u[7] + σ[5]*u[5]*u[8] + σ[6]*u[1]*u[3] + σ[7]*u[3]*u[9];
    du[3] = -u[3]*k[3]/Re + σ[8]*(u[4]*u[7]+u[5]*u[6]) + σ[9]*u[4]*u[8]; 
    du[4] = -u[4]*k[4]/Re + σ[10]*u[1]*u[5] + σ[11]*u[2]*u[6] + σ[12]*u[3]*u[7] + σ[13]*u[3]*u[8] + σ[14]*u[5]*u[9]; 
    du[5] = -u[5]*k[5]/Re + σ[15]*u[1]*u[4] + σ[16]*u[2]*u[7] + σ[17]*u[2]*u[8] + σ[18]*u[4]*u[9] + σ[19]*u[3]*u[6]; 
    du[6] = -u[6]*k[6]/Re + σ[20]*u[1]*u[7] + σ[21]*u[1]*u[8] + σ[22]*u[2]*u[4]+ σ[23]*u[3]*u[5] + σ[24]*u[7]*u[9] + σ[25]*u[8]*u[9]
    du[7] = -u[7]*k[7]/Re + σ[26]*(u[1]*u[6]+u[6]*u[9]) + σ[27]*u[2]*u[5] + σ[28]*u[3]*u[4]
    du[8] = -u[8]*k[8]/Re + σ[29]*u[2]*u[5] + σ[30]*u[3]*u[4] 
    du[9] = -u[9]*k[9]/Re + σ[31]*u[2]*u[3] + σ[32]*u[6]*u[8] 
end

function compute_E9D(di::Dict)
    @unpack res, Re = di
    p = E9DParameters(; Re = Re)
    ds = ContinuousDynamicalSystem(E9D!, zeros(9), p, (J,z0, p, n) -> nothing)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    yg = range(-5, 5; length = 10001)
    grid = ntuple(x -> yg, 9)
    mapper = AttractorsViaRecurrences(ds, grid; sparse = true, Δt = .01,   
        mx_chk_fnd_att = 300, 
        mx_chk_loc_att = 100, safety_counter_max = Int(1e7), diffeq)
    u0(x,y) = [x, 0.0511, -0.0391, 0.0016, y, 0.126, 0., 0., 0.]
    y1r = range(-1, 1, length = res)
    y2r = range(-1, 1, length = res)
    ics = [ u0(y1,y2) for y1 in y1r, y2 in y2r]
    bsn = [ mapper(u) for u in ics]
    grid = (y1r,y2r)
    return @strdict(bsn, grid, Re, res)
end


function print_basins(w,h,cmap, Re, res)
    params = @strdict res Re
    @time data, file = produce_or_load(
        datadir("basins"), params, compute_E9D;
        prefix = "eckhardt", storepatch = false, suffix = "jld2", force = true
    )
    @unpack bsn, grid = data
    xg, yg = grid
    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"$y$", xlabel = L"x", yticklabelsize = 30, 
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
    save(string(projectdir(), "/plots/eckhardt_",res,".png"),fig)
end




print_basins(600,600, nothing, 407., 300)
# f,a,r = continuation_E9D()
# plot_filled_curves(f,r,"tst.png")
