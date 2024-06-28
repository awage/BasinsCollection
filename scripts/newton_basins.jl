using DrWatson
@quickactivate 
using Attractors
using CairoMakie
using LaTeXStrings
using ForwardDiff: derivative
using LinearAlgebra:norm

func_list = [x -> (x*x - 1) * (x*x + 1),
x -> x*x*x - 1,
x -> x^12 - 1,
x -> (x*x - 4)*(x + 1.5)*(x - 0.5),
x -> (x+2)*(x+1.5)^2*(x-0.5)*(x-2),
x -> sin(x),
x -> (x - 1)^3 + 4 * (x-1)^2 - 10,
x -> sin(x-14/10)^2 - (x - 14/10)^2 + 1,
x -> x*x - exp(x) - 3x + 2,
x -> cos(x-3/4) - x + 3/4,
x -> (x + 1)^3 - 1,
x -> (x-2)^3 - 10,
x -> (x + 5/4) * exp((x + 5/4)*(x + 5/4)) - sin((x + 5/4))^2 + 3 * cos((x + 5/4)) + 5,
x -> (x + sin(2/x) * x*x)]

string_list = [L"f(x) = (x^2 - 1)(x^2 + 1)",
L"f(x) = x^3 - 1",
L"f(x) = x^{12} - 1",
L"f(x) = (x^2 - 4)(x + 1.5)(x - 0.5)",
L"f(x) =(x+2)(x+1.5)^2 (x-0.5)(x-2)",
L"f(x) = \sin(x)",
L"f(x) = (x - 1)^3 + 4(x-1)^2 - 10",
L"f(x) = \sin(x-14/10)^2 - (x - 14/10)^2 + 1",
L"f(x) = x^2 - e^x - 3x + 2",
L"f(x) = \cos(x-3/4) - x + 3/4",
L"f(x) = (x + 1)^3 - 1",
L"f(x) = (x-2)^3 - 10",
L"f(x) = (x + 5/4)~  e^{(x + 5/4)^2} - \sin(x + 5/4)^2 + 3\cos(x + 5/4) + 5",
L"f(x) = (x + sin(2/x)  x^2)"]


function ∂f(f)
# Warning. This trick works only for holomorphic functions. 
    ∂f∂z(z) = (g(x) = f(x+z); derivative(x->real(g(x)),0) 
            + im * derivative(x->imag(g(x)),0))
    return ∂f∂z
end

function N_map(z, f, ∂f∂z)
    dz = f(z)/∂f∂z(z)
    return  z - dz
end

function beta_map(f)
    ∂f∂z = ∂f(f)
    N(z) = N_map(z, f, ∂f∂z)
    function N_β(z1, p, n)
        β = p[1]
        z = z1[1] + im * z1[2]
        N_z = N(z)
        z_new =  N_z - β * f(N_z)/∂f∂z(z)
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
    @unpack i, res = di
    ε = 1e-14; max_it = 50; β = -1
    N_β = beta_map(func_list[i])
    ds = DiscreteDynamicalSystem(N_β, [0.1, 0.2], [β])
    xg = yg = range(-10, 10; length = 10000)
    grid = (xg, yg)
    # We set up a mapper so that we can identify roots automatically  
    mapper_beta = AttractorsViaRecurrences(ds, (xg, yg);
            sparse = true, consecutive_recurrences = 3000
    )
    xg = yg = range(-2, 2; length = res)
    grid = (xg, yg)
    bsn = zeros(Int8, res,res); 
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

# cmap = ColorScheme([RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )
cmap = :Spectral_7

for i in 1:14
    res = 1000
    params = @strdict  res i
    print_fig(params, string("newton_",i),compute_basins_newton; ylab = L"\Im{(z)}", xlab = L"\Re{(z)}", cmap)
# print_fig(600, 600, cmap, 3, 2000) 
end
