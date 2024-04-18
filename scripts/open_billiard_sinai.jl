using DrWatson 
@quickactivate
using DynamicalBilliards
using LaTeXStrings
using CairoMakie

function escapewall!(p::AbstractParticle{T}, bd::Billiard{T}, t)::T where {T<:AbstractFloat}

    ei = DynamicalBilliards.escapeind(bd)

    totalt = zero(T); count = zero(t); t_to_write = zero(T)
    i = 0
    while count < t

        i, tmin, pos, vel = bounce!(p, bd)
        t_to_write += tmin

        if DynamicalBilliards.isperiodic(i, bd)
            continue
        else
            totalt += t_to_write
            i ∈ ei &&  break # the collision happens with a Door!

            # set counter
            count += DynamicalBilliards.increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end#time, or collision number, loop
    if count ≥ t
        warning && warn("Particle did not escape within max-iter window.")
        return T(Inf)
    end
    return i
end


#PHYSICAL REVIEW A VOLUME 38, NUMBER 2 JULY 15, 1988
#Fractal boundaries for exit in Hamiltonian dynamics
#Siegfried Bleher, ' Celso Grebogi, Edward Ott, and Reggie Brow
function get_billiard(Δ)
    bd = Obstacle{Float64}[]
    a = 0.1; b = 0.2; L = 4.; r = 1; 
    if Δ ≥ L/2
        error("Δ too big")
    end
    
    a = (L - 2Δ)/4
    b = 2*a 

    wi = FiniteWall([-L/2,-L/2], [-L/2,L/2], [1.,0], "muro izdo")
    wd = FiniteWall([L/2,-L/2], [L/2,L/2], [-1.,0], "muro dcho")
    wu = FiniteWall([-L/2,L/2], [L/2,L/2], [0,-1], "muro arriba")

    v = [a Δ b Δ a]
    vs = cumsum(vec(v))
    wd1 = FiniteWall([-L/2,-L/2], [-L/2+vs[1], -L/2], [0, 1.], "muro b1")
    wd2 = FiniteWall([-L/2+vs[1],-L/2], [-L/2+vs[2], -L/2], [0, 1.], true, "muro b2 abierto")
    wd3 = FiniteWall([-L/2+vs[2],-L/2], [-L/2+vs[3], -L/2], [0, 1.], "muro b3")
    wd4 = FiniteWall([-L/2+vs[3],-L/2], [-L/2+vs[4], -L/2], [0, 1.], true, "muro b4 abierto")
    wd5 = FiniteWall([-L/2+vs[4],-L/2], [-L/2+vs[5], -L/2], [0, 1.], "muro b5")
    d1 = Disk([0,0], r, "disco central")
    push!(bd, wi, wd, wu, wd1, wd2, wd3, wd4, wd5, d1)
    return Billiard(bd) 
end

function compute_exit_open_sinai(di::Dict)
    @unpack Δ, res = di
    br = get_billiard(Δ)
    # br = Billiard(bd)
    x0 = range(-1, 1, length = res) 
    θ0 = range(-pi+pi/2,pi+pi/2,length = res)
    u(x,θ) = Particle([x, 1.25, θ])
    bsn = [ escapewall!(u(x,θ), br, 10999) for x in x0, θ in θ0] 
    θ0 = range(-pi,pi,length = res)
    grid = (x0,θ0)
    return @strdict(bsn, grid, Δ, res)
end

Δ = 0.8; res = 700
params = @strdict Δ res
print_fig(params, "open_sinai", compute_exit_open_sinai)
