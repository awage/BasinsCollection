using DrWatson
@quickactivate
using Attractors
using CairoMakie
using LaTeXStrings
# International Journal of Bifurcation and ChaosVol. 04, No. 02, pp. 343-381 (1994) 
# BASIN BIFURCATIONS OF TWO-DIMENSIONAL NONINVERTIBLE MAPS: FRACTALIZATION OF BASINS
# C. MIRA , D. FOURNIER-PRUNARET, L. GARDINI , H. KAWAKAMI, and J.C. CATHALA
# https://doi.org/10.1142/S0218127494000241
function quadratic_map!(dz, z, p, n)
    xn, yn = z
    a, b = p
    dz[1] = a*xn + yn
    dz[2] = xn^2 + b 
    return
end


function compute_mira(di)
    @unpack a, b, res= di
    u0 = [0., 0.6]
    ds = DeterministicIteratedMap(quadratic_map!, [1.0, 0.0], [a, b])
    xgg = range(-10, 10, length = 40001)
    ygg = range(-10, 10, length = 40001)
    grid = (xgg, ygg)
    mapper = AttractorsViaRecurrences(ds, grid)
            # mx_chk_fnd_att = 3000,
            # mx_chk_loc_att = 3000, sparse = true)
    xg = range(-2.05, 2.05, length = res)
    yg = range(-3, 3, length = res)
    grid = (xg,yg)
    bsn, att = basins_of_attraction(mapper, grid; show_progress = true)
    return @strdict(bsn, att,grid)
end

a = -0.42
b = -1.32
params = @strdict res a b
cmap = ColorScheme([RGB(1,1,1),  RGB(0.9,0.2,0.1)] )
print_fig(params, "mira_map", compute_mira; xlab = L"x", ylab = L"y", force = false, cmap)
