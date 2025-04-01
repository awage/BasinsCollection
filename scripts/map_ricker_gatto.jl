using DrWatson
@quickactivate
using Attractors
using CairoMakie
using LaTeXStrings
using Colors,ColorSchemes
include(srcdir("print_fig.jl"))
# Chaos, Solitons and Fractals 12 (2001) 301±311
# Dynamics with riddled basins of attraction in models of
# interacting populations
# Bernard Cazelles
function ricker_gatto!(dz, z, p, n)
    xn, yn = z
    r1 = 22; s1 = 0.007815; r2 = 2.95; s2 = 0.5; 
    dz[1] = xn*(r1*exp(-xn - yn) + s1)
    dz[2] = yn*(r2*exp(-xn - yn) + s2)
    return
end


function compute_RG(di)
    @unpack res= di
    ds = DeterministicIteratedMap(ricker_gatto!, rand(2))
    xgg = range(0, 10, length = 100001)
    ygg = range(0, 1, length = 100001)
    mapper = AttractorsViaRecurrences(ds, (xgg,ygg); 
                    consecutive_recurrences = 1000, 
                    consecutive_attractor_steps = 10
                                     )
    xg = range(3, 6, length = res)
    yg = range(0, 0.1, length = res)
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att,grid)
end

res = 2000
params = @strdict res 
cmap = ColorScheme([RGB(0.1,0.1,0.1), RGB(1,1,1)] )
print_fig(params, "ricker_gatto", compute_RG; xlab = L"x", ylab = L"y", force = false, cmap)
