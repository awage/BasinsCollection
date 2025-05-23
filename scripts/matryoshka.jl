using DrWatson
@quickactivate
using OrdinaryDiffEqVerner
using Attractors
using CairoMakie
using ProgressMeter
using LaTeXStrings
include(srcdir("print_fig.jl"))

# Matryoshka multistability: Coexistence of an infinite number of exactly
# self-similar nested attractors in a fractal phase space
# Artur Karimov, Ivan Babkin, Vyacheslav Rybin, Denis Butusov 
# Chaos, Solitons and Fractals 187 (2024) 115412
# https://doi.org/10.1016/j.chaos.2024.115412
function Fr(y) 
    d = 1.; m = 1.; 
    P = 1.23; R = 2.
    ay = abs(y)
    while true
        if ay < d 
            d = d/R
            m = m/R
        elseif ay > 2*d
            d = R*d 
            m = R*m
        else 
            break
        end
    end
    ε = 0.01
    if d > ε
        if ay < P*d 
            b = -m*(R-P*R+1)/(R*(P-1))
            k = -m/(R*d*(1-P))
        else
            b = -m*(R-P*R+1)/(P-R)
            k = m*(-R^2+R+1)/(R*d*(P-R))
        end
        return k*ay + b
    else
        return 0 
    end
end

function matryoshka!(du, u, p, t)
    x,y,z = u
    a = 1.9; b = -1.8; c = 3.9
    du[1] = a*z
    du[2] = b*y + z
    du[3] = -x + y + c*Fr(y)
    return nothing
end

function compute_matryoshka(di::Dict)
    @unpack  res = di
    diffeq = (reltol = 1e-6,  alg = Vern9())
    df = CoupledODEs(matryoshka!,rand(3), [])
    yg = range(-1000, 1000; length = 8001)
    grid = ntuple(x -> yg, 3)
    pow = 7; d = 1000
    # using an irregular grid because of the accumumlation of 
    # attractors near the origin. 
    xg = range(-(100^(1/pow)), 1000^(1/pow); length = d).^pow
    yg = range(-(200^(1/pow)), 200^(1/pow); length = d).^pow
    zg = range(-(200^(1/pow)), 200^(1/pow); length = d).^pow
    grid = (xg, yg, zg)
    mapper = AttractorsViaRecurrences(df, grid; Δt = 0.01, Ttr = 10)
    xg = range(-5,80,length = res)
    yg = range(-20,30,length = res)
    bsn = @showprogress [ mapper([x, 0., y]) for x in xg, y in yg] 
    grid = (xg, yg)
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid,  res)
end


res = 1200
params = @strdict res
print_fig(params, "matryoshka", compute_matryoshka; force = false, ylab = L"z")

# data, file = produce_or_load(
#     datadir("basins"), params, compute_matryoshka;
#     prefix = "matryoshka", storepatch = false, suffix = "jld2", force = true
# )
# @unpack bsn, grid,att = data

