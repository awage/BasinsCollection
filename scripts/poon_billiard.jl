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
        @warn("Particle did not escape within max-iter window.")
        return T(Inf)
    end
    return i
end

function get_billiard_poon(w)
    bd = Obstacle{Float64}[]
    wd1 = FiniteWall([-0.5,0], [0, sqrt(3)/2], [cos(-pi/6), sin(-pi/6)], true, "muro b1")
    wd2 = FiniteWall([0, sqrt(3)/2], [0.5, 0], [cos(pi+pi/6), sin(pi+pi/6)], true, "muro b2 abierto")
    wd3 = FiniteWall([-0.5, 0], [0.5, 0], [0,1.], true, "muro b3")
    r = (1 - w)/2 
    d1 = Disk([-0.5,0], r, "disco izqdo")
    d2 = Disk([0.5,0], r, "disco dcho")
    d3 = Disk([0,sqrt(3)/2], r, "disco arrib")

    we1 = FiniteWall([-3.,-3.], [-3., 3.], [1., 0.], true, "muro b1")
    we2 = FiniteWall([-3.,3.], [3., 3.], [0,-1], true, "muro b1")
    we3 = FiniteWall([3.,3.], [3., -3.], [-1, 0], true, "muro b1")
    we4 = FiniteWall([3.,-3.], [-3., -3.], [0, 1], true, "muro b1")

    push!(bd, wd1, wd2, wd3, d1, d2, d3, we1, we2, we3, we4)

    br = Billiard(bd)
    return br
end

function compute_exit_poon_disks(di::Dict)
    @unpack ww, res = di
    br = get_billiard_poon(ww)
    y0 = range(-0.2, 0.2, length = res) 
    θ0 = range(1.4, 2,length = res)
    u(y,θ) = Particle([y, -0.5,  θ])
    bsn = [ escapewall!(u(y,θ), br, 20000000) for y in y0, θ in θ0] 
    bsn[bsn .>= 7] .= 0
    grid = (y0,θ0)
    return @strdict(bsn, grid, ww, res)
end


ww = 0.1; res = 1400
params = @strdict ww res
print_fig(params, "poon_disks", compute_exit_poon_disks)
