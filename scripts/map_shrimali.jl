using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using ProgressMeter
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))

  # title={Basin bifurcations in quasiperiodically forced coupled systems},
  # author={Shrimali, Manish Dev and Prasad, Awadhesh and Ramaswamy, Ram and Feudel, Ulrike},
  # journal={Physical Review E—Statistical, Nonlinear, and Soft Matter Physics},
  # volume={72},
  # number={3},
  # pages={036215},
  # year={2005},
function map_shrimali!(dz, z, p, n)
x,y,θ = z
α, β, ϵp = p
τ = (sqrt(5)-1)/2
ϵ = ϵp*(4/α -1)  
dx = α*(1 + ϵ*cos(2π*θ))*x*(1 - x) + β*(y - x)
dy = α*(1 + ϵ*cos(2π*θ))y*(1 - y) + β*(x - y)
dθ = mod(θ + τ,1) 
dz[1] = dx; dz[2] = dy; dz[3] = dθ
    return
end

function compute_shrimali(di::Dict)
    @unpack α, ϵp, β, res = di
    ds = DeterministicIteratedMap(map_shrimali!, rand(3), [α, β, ϵp])
    yg = xg =  range(0, 1., length = 5001)
    projection = [1,2]; complete_state = [0.]
    pinteg = ProjectedDynamicalSystem(ds, projection, complete_state)
    mapper = AttractorsViaRecurrences(pinteg, (xg, yg))
    yg = xg = range(0., 1., length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att, grid, res)
end

res = 1200; α = 3.2511; ϵp = 0.5  ; β = 0.01 
params = @strdict res α β ϵp 
cmap = ColorScheme([RGB(1,1,1),  RGB(1,1,1),  RGB(0.9,0.2,0.1),  RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
print_fig(params, "shrimali", compute_shrimali; xlab = L"x", ylab = L"y", force = false, cmap)
# att =  get_att(params, "nash_equilibrium", compute_nash; force = true)

script_name = basename(Base.source_path())

println("Current script name: $script_name")
