using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter


function compute_thomas(di::Dict)
    @unpack b, res = di
    ds = Systems.thomas_cyclical(b = b)
    xg = yg = zg = range(-7, 7,length=10000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, maximum_iterations = Int(1e7), show_progress = true)
    y1 = y2 = range(-5, 5, length = res)
    bsn = @showprogress [ mapper([x,y,0.]) for x in y1, y in y2]
    att = mapper.bsn_nfo.attractors
    grid = (y1,y2)
    return @strdict(bsn, att, grid, b, res)
end


b=0.1665; #res = 500
params = @strdict res b
print_fig(params, "thomas", compute_thomas; ylab = L"\dot{x}", xlab = L"x")
