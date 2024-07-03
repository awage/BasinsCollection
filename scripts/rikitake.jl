using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter
# This system is problematic: very very long transient (t > 2000 sometimes) that can be mistaken with attractors.




function compute_rikitake(di::Dict)
    @unpack μ, α, res = di
    ds = Systems.rikitake(μ = μ, α = α)
    xg = yg = zg = range(-5,5,length=10000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg); sparse = true,    
        mx_chk_fnd_att = 10000,
        mx_chk_loc_att = 10000, maximum_iterations = Int(1e7), show_progress = true)
    y1 = y2 = range(-2.5, 2.5, length = res)
    bsn = @showprogress [ mapper([x,y,0.]) for x in y1, y in y2]
    att = mapper.bsn_nfo.attractors
    grid = (y1,y2)
    return @strdict(bsn, att, grid, μ, α, res)
end



μ = 0.47; α = 1.; #res = 700
params = @strdict res μ α
print_fig(params, "rikitake", compute_rikitake; ylab = L"\dot{x}", xlab = L"x")
