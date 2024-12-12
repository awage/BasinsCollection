using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie

# C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983
function grebogi_map(dz,z, p, n)
    θ = z[1]; x = z[2]
    J₀=0.3; a=1.32; b=0.9;
    dz[1]= θ + a*sin(2*θ) - b*sin(4*θ) -x*sin(θ)
    dz[1] = mod(dz[1],2π) # to avoid problems with attracto at θ=π
    dz[2]=-J₀*cos(θ)
    return
end


function compute_grebogi(di::Dict)
    @unpack res = di
    ds = DeterministicIteratedMap(grebogi_map,[1., -1.], [] )
    θ = range(0, 2π, length = 2000)
    xg = range(-0.5, 0.5, length = 2000)
    mapper = AttractorsViaRecurrences(ds, (θ,xg))
    θ = range(0, 2π, length = res)
    xg = range(-0.5, 0.5, length = res)
    bsn, att = basins_of_attraction(mapper, (θ,xg); show_progress = true)
    grid = (θ,xg)
    return @strdict(bsn, att, grid, b, res)
end


params = @strdict res 
cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
print_fig(params, "grebogi", compute_grebogi; ylab = L"x", xlab = L"\theta", cmap)
