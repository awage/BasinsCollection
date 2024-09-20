using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
using OrdinaryDiffEq:Vern9
using ProgressMeter

"""
    thomas_cyclical(u0 = [1.0, 0, 0]; b = 0.2)
```math
\\begin{aligned}
\\dot{x} &= \\sin(y) - bx\\\\
\\dot{y} &= \\sin(z) - by\\\\
\\dot{z} &= \\sin(x) - bz
\\end{aligned}
```
Thomas' cyclically symmetric attractor is a 3D strange attractor originally proposed
by René Thomas[^Thomas1999]. It has a simple form which is cyclically symmetric in the
x,y, and z variables and can be viewed as the trajectory of a frictionally dampened
particle moving in a 3D lattice of forces.
For more see the [Wikipedia page](https://en.wikipedia.org/wiki/Thomas%27_cyclically_symmetric_attractor).

Reduces to the labyrinth system for `b=0`, see
See discussion in Section 4.4.3 of "Elegant Chaos" by J. C. Sprott.

[^Thomas1999]:
    Thomas, R. (1999). *International Journal of Bifurcation and Chaos*,
    *9*(10), 1889-1905.
"""
function thomas_rule(u, p, t)
    x,y,z = u
    b = p[1]
    xdot = sin(y) - b*x
    ydot = sin(z) - b*y
    zdot = sin(x) - b*z
    return SVector{3}(xdot, ydot, zdot)
end

function compute_thomas(di::Dict)
    @unpack b, res = di
    ds =  CoupledODEs(thomas_rule, rand(3), [b])
    xg = yg = zg = range(-7, 7,length=10000)
    mapper = AttractorsViaRecurrences(ds, (xg,yg,zg))
    y1 = y2 = range(-5, 5, length = res)
    bsn = @showprogress [ mapper([x,y,0.]) for x in y1, y in y2]
    att = mapper.bsn_nfo.attractors
    grid = (y1,y2)
    return @strdict(bsn, att, grid, b, res)
end


b=0.1665; #res = 500
params = @strdict res b
print_fig(params, "thomas", compute_thomas; ylab = L"y", xlab = L"x")
