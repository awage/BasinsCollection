using DrWatson
@quickactivate
using Attractors
using LaTeXStrings
using CairoMakie
include("../src/print_fig.jl")

function Ms(x)
    y=0.
     if x ≤ 0; y=9*x/(4-5*x); end
     if 0 < x ≤ 4/9; y=9/4*x; end
     if 4/9 < x ≤ 5/9; y=81/4*(x-x^2)-4; end
     if  x > 5/9; y=9/4*(1-x); end
     return y
 end

function chaotic_map(dz,z, p, n)
    xn = z[1]; yn = z[2]
    λ=p[1];
    dz[1]=Ms(xn)
    dz[2]=λ*yn+sin(2*π*xn)
    return
end


function escape_function(y)
    if y > 10000
        return 1
    elseif y < -10000
        return 2
    end
    return 0
end

function  get_color(integ,x,y)
    reinit!(integ,[x,y])
    while escape_function(integ.u[2]) == 0
        step!(integ)
        integ.t > 2000 && break
    end
    return escape_function(integ.u[2]);
end

function compute_sporadical(di)
    @unpack res = di
    xg=range(0.,1.,length=600)
    yg=range(-2.,6.,length=600)
    ds = DeterministicIteratedMap(chaotic_map,[1., 0.], [1.1])
    integ  = integrator(ds)
    integ.p[1] = 1.1
    bsn = [get_color(integ,x,y) for x in xg, y in yg]
    grid = (xg,yg)
    return @strdict(bsn, grid)
end



# res = 500
params = @strdict res
print_fig(params, "sporadical", compute_sporadical; xlab = L"x", ylab = L"y")


