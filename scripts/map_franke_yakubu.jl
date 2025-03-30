using DrWatson
@quickactivate
using Attractors
using CairoMakie
using LaTeXStrings

# Chaos, Solitons and Fractals 12 (2001) 301±311
# Dynamics with riddled basins of attraction in models of
# interacting populations
# Bernard Cazelles
function Franke_Yakubu!(dz, z, p, n)
    xn, yn = z
    r1 = 2.825; s1 = 20.25; c1 = 1.2; c2 = 0.1
    dz[1] = xn*exp(r1 - s1*(xn + yn))
    dz[2] = c1*yn/(c2 + xn + yn) 
    return
end


function compute_FY(di)
    @unpack res= di
    ds = DeterministicIteratedMap(Franke_Yakubu!, rand(2))
    xgg = range(0, 100, length = 40001)
    ygg = range(0, 2, length = 40001)
    mapper = AttractorsViaRecurrences(ds, (xgg,ygg))
    xg = range(0, 30, length = res)
    yg = range(0, 2, length = res)
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att,grid)
end

res = 1200
params = @strdict res 
# cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
print_fig(params, "franke_yakubu", compute_FY; xlab = L"x", ylab = L"y", force = true, cmap)
