using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter


# International Journal of Bifurcation and Chaos, Vol. 26, No. 14 (2016) 1650233 (11 pages)
#DOI: 10.1142/S0218127416502333
#Crisis in Amplitude Control Hides in Multistability
#Int. J. Bifurcation Chaos 2016.26.
# Chunbiao Li, Julien Clinton Sprott, Hongyan Xing. 
function li_sprott!(du, u, p, t)
    @inbounds begin
    a = p[1]; b = p[2];
    x, y, z = u
    du[1] = y + y*z
    du[2] = y*z - a*x*z
    du[3] = b*z^2 - y^2
    end
end


function compute_li_sprott(di::Dict)
    @unpack a, b, res = di
    ds = CoupledODEs(li_sprott!, rand(3), [a,b])
    xg = yg = zg = range(-5,5,length=10000)
    diffeq = (alg = Vern9(), reltol = 1e-9, maxiters = 1e8)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); sparse = true,    
        consecutive_recurrences = 1000,
        attractor_locate_steps = 1000, maximum_iterations = Int(1e8), show_progress = true)
    x1 = range(-1, 1, length = res) 
    y1 = range(-5, 5, length = res)
    bsn = @showprogress [ mapper([x,y,-1.]) for x in x1, y in y1]
    att = mapper.bsn_nfo.attractors
    # bsn, att = basins_of_attraction(mapper, (y1,y2); show_progress = true)
    grid = (x1,y1)
    return @strdict(bsn, att, grid, res)
end

a = 13.; b = 0.55; #res = 500
params = @strdict res a b
cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(0.7,0.7,0.7), RGB(1,0,0)] )
print_fig(params, "li_sprott", compute_li_sprott; cmap) 


