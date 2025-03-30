using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
include(srcdir("print_fig.jl"))

# James C. Alexander, Brian R. Hunt, Ittai Kan and James A. Yorke
# Ergodic Theory and Dynamical Systems / Volume 16 / Issue 04 / August 1996, pp 651 - 662
# DOI: 10.1017/S0143385700009020,
function F(x,y,λ)     
    z = x + im*y
    dz = (z^2 - ( 1+ im*λ)*conj(z))
    dz = (dz^2 - ( 1+ im*λ)*conj(dz))
    return real(dz), imag(dz)
end

function triangular_map(dz, z, p, n)
    xn = z[1]; yn = z[2]
    λ = p
    dz[1], dz[2] = F(xn,yn,λ)
    return
end

function compute_triangular(di::Dict)
    @unpack  λ, res = di
    ds = DeterministicIteratedMap(triangular_map, [1.0, 0.0], λ)
    yg = xg = range(-1.5, 1.5, length = 105000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg);     
        Ttr = 100,
        consecutive_recurrences = 1000,
        attractor_locate_steps = 1000,
       consecutive_attractor_steps = 5)
    yg = xg = range(-2, 2, length = res)
    # For this system we have to seed the attractors since there are some crossing 
    # regions. See the paper.
    mapper([0.6758665951691236, 0.14192158197150856])
    mapper( [1.2082683535288237, 0.51435685])
    mapper( [-0.8864174534443858, 0.5066063660814969])
    mapper( [-0.18423411840724807, -0.709610846430985])
    mapper( [0.3938998475517044, -0.34645915094681234])
    mapper( [-0.48700615019363475, 0.5143568499999996])
    bsn, att = basins_of_attraction(mapper, (xg,yg); show_progress = true)
    grid = (xg, yg)
    return @strdict(bsn, att, grid, res)
end

res = 3000
λ = 1.0287137
params = @strdict res λ
# cmap = ColorScheme([RGB(1,1,1), RGB(0,1,0), RGB(0.34,0.34,1), RGB(1,0.46,0.46), RGB(0.1,0.1,0.1) ] )
print_fig(params, "triangular_img", compute_triangular; force = true) 
att = get_att(params, "triangular_img", compute_triangular) 

