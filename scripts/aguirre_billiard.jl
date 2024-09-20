using DrWatson 
@quickactivate
using DynamicalBilliards
using LaTeXStrings
using CairoMakie
using ProgressMeter

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

# Limit of small exits in open Hamiltonian systems
# Jacobo Aguirre and Miguel A. F. Sanjuan
# Phys Rev E 67 2003 ́
function get_billiard(w)
    bd = Obstacle{Float64}[]
    wd1 = FiniteWall([-0.5,0], [0, sqrt(3)/2], [cos(-pi/6), sin(-pi/6)], true, "muro b1")
    wd2 = FiniteWall([0, sqrt(3)/2], [0.5, 0], [cos(pi+pi/6), sin(pi+pi/6)], true, "muro b2 abierto")
    wd3 = FiniteWall([-0.5, 0], [0.5, 0], [0,1.], true, "muro b3")
    r = (1 - w)/2 
    d1 = Disk([-0.5,0], r, "disco izqdo")
    d2 = Disk([0.5,0], r, "disco dcho")
    d3 = Disk([0,sqrt(3)/2], r, "disco arrib")
    push!(bd, wd1, wd2, wd3, d1, d2, d3)

    br = Billiard(bd)
    return br
end

function compute_exit_three_disk(di::Dict)
    @unpack ww, res = di
    br = get_billiard(ww)
    x0 = range(-ww/2, ww/2, length = res) 
    θ0 = range(0,pi,length = res)
    u(x,θ) = Particle([x, 0.0001, θ])
    bsn = @showprogress [ escapewall!(u(x,θ), br, 20000000) for x in x0, θ in θ0] 
    grid = (x0,θ0)
    return @strdict(bsn, grid, ww, res)
end


ww = 0.001; #res = 400
params = @strdict ww res
print_fig(params, "open_disks", compute_exit_three_disk; xlab = L"x_0", ylab = L"\theta_0")
