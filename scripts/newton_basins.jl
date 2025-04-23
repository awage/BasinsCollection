using DrWatson
@quickactivate 
using Attractors
using CairoMakie
using ProgressMeter
using Colors; using ColorSchemes
using LaTeXStrings
using ForwardDiff: derivative
using LinearAlgebra:norm
include(srcdir("print_fig.jl"))


function ∂f(f)
# Warning. This trick works only for holomorphic functions. 
    ∂f∂z(z) = (g(x) = f(x+z); derivative(x->real(g(x)),0) 
            + im * derivative(x->imag(g(x)),0))
    return ∂f∂z
end


function Newton_map(f)
    ∂f∂z = ∂f(f)
    function N_β(z1, p, n)
        z = z1[1] + im * z1[2]
        z_new = z - f(z)/∂f∂z(z)
        return SVector(real(z_new), imag(z_new))
    end
    return N_β
end

function _get_iterations!(ds, ε, max_it)
    xn_1 = get_state(ds) 
    step!(ds)
    xn = get_state(ds) 
    k = 1
    # stopping criterion is ∥x_n - x_{n-1}∥ ≤ ε
    while norm(xn - xn_1) > ε
        (k > max_it) && break 
        xn_1 = xn
        step!(ds)
        xn = get_state(ds) 
        k += 1
    end
    return k
end


function compute_basins_newton(di::Dict)
    @unpack f, res = di
    ε = 1e-14; max_it = 50; β = 0.
    N = Newton_map(f)
    ds = DiscreteDynamicalSystem(N, [0.1, 0.2])
    xg = yg = range(-10, 10; length = 10000)
    grid = (xg, yg)
    # We set up a mapper so that we can identify roots automatically  
    mapper_beta = AttractorsViaRecurrences(ds, (xg, yg);
            sparse = true, consecutive_recurrences = 3000
    )
    xg = yg = range(-2, 2; length = res)
    grid = (xg, yg)
    bsn = zeros(Int, res,res); 
@showprogress for (i,x) in enumerate(xg), (j,y) in enumerate(yg) 
        set_state!(ds, [x,y])
        n = _get_iterations!(ds,ε,max_it)
        if n > max_it
            # the alg. did not converge
            bsn[i,j] = -1
        else
            # We identify the root with the mapper.
            bsn[i,j] = mapper_beta([x,y])
        end
    end
    att = extract_attractors(mapper_beta)
    return @strdict(bsn, att, grid,  res)
end

# cmap = ColorScheme([RGB(0,0,0), RGB(1,0,0), RGB(0,1,0), RGB(0,0,1), RGB(1,1,0), RGB(1,0,1), RGB(0,1,1)] )
cmap = :flag
f(x) = (x + sin(2/x) * x*x)
params = @strdict  res f
print_fig(params,"newton_",compute_basins_newton; ylab = L"\Im{(z)}", xlab = L"\Re{(z)}", cmap, force = false)
