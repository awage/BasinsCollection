using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes


# (M+m)*ddxs + α*cos(ϕ1)*ddϕ1 + m*l*cos(ϕ2)*ddϕ2 = α*dϕ1^2*sin(ϕ1) + m*l*dϕ2^2*sin(ϕ2) - kx*xs - dx*dxs

# -β*sin(ϕ1)*ddϕ1 - m*l*sin(ϕ2)*ddϕ2 + (M+m)*ddys = β*cos(ϕ1)*dϕ1^2 + m*l*cos(ϕ2)*dϕ2^2 + g*(m + M) - ky*ys - dy*dys

# m*l*cos(ϕ2)*ddx_s + m*l*lc*cos(ϕ1-ϕ2)*ddϕ1 + BC*ddϕ2 - m*l*sin(ϕ2)*ddy_s = m*l*lc*dϕ1^2*sin(ϕ1-ϕ2) - m*g*l*sin(ϕ2) - dc*(dϕ2 - dϕ1)

# α*cos(ϕ1)*ddxs + γ*ddϕ1 + m*l*lc*cos(ϕ1-ϕ2)*ddϕ2 - β*sin(ϕ1)*ddys = -m*l*lc*sin(ϕ1-ϕ2)*dϕ2^2 - α*g*sin(ϕ1) - db*dϕ1 - dc*(dϕ1 - dϕ2) - dbf*2/pi*atan(10^5*dϕ1) - dcf*2/pi*atan(10^5*(dϕ1 - dϕ2)) + Mt

function bell_yoke!(du, u, p, t)
M = 147.55            
BB = 23.86             
m = 5.49              
BC = 0.52              
L  = 0.17              
l   = 0.27              
lc = 0.13
kx = 813770.00        
ky = 10939000.00     
dx = 260.30            
dy = 954.00            
g = 9.81              
db = 0.22              
dc = 0.10              
dbf = 0.10              
dcf = 0.25              
r  = 0.255             
k = 0.125             
α = 0.497             
β = 0.424             
α = M*L + M*lc
β = m*lc + M*L
γ = BB + m*lc^2
A = 15.; Tmax = 20.; ω = 7.5 
Mt(ϕ,dϕ) = abs(ϕ) ≤ π/A ? Tmax*sign(dϕ)*cos(ω*ϕ) : 0.

ϕ1, ϕ2, xs, ys, dϕ1, dϕ2, dxs, dys  = u
du[1] = dϕ1; du[2] = dϕ2; du[3] = dxs; du[4] = dys

AA = [
α*cos(ϕ1) m*l*cos(ϕ2) (M+m) 0; 
-β*sin(ϕ1) -m*l*sin(ϕ2) 0 M+m; 
m*l*lc*cos(ϕ1-ϕ2) BC m*l*cos(ϕ2) -m*l*sin(ϕ2); 
γ m*l*lc*cos(ϕ1-ϕ2) α*cos(ϕ1) -β*sin(ϕ1)]


b1 = α*dϕ1^2*sin(ϕ1) + m*l*dϕ2^2*sin(ϕ2) - kx*xs - dx*dxs; 
b2 = β*cos(ϕ1)*dϕ1^2 + m*l*cos(ϕ2)*dϕ2^2 + g*(m + M) - ky*ys  -  dy*dys;
b3 = m*l*lc*dϕ1^2*sin(ϕ1 - ϕ2) - m*g*l*sin(ϕ2) - dc*(dϕ2 - dϕ1);
b4 = -m*l*lc*sin(ϕ1 - ϕ2)*dϕ2^2 - β*g*sin(ϕ1) - db*dϕ1 - dc*(dϕ1 - dϕ2) - dbf*2/pi*atan(10^5*dϕ1) - dcf*2/pi*atan(10^5*(dϕ1 - dϕ2))+Mt(ϕ1,dϕ1)

du[5:8] = AA\[b1;b2;b3;b4]

end


function compute_bell_yoke(di::Dict)
    @unpack res = di
    diffeq = (reltol = 1e-6,  alg = Vern9(), maxiter = 1e6)
    df = CoupledODEs(bell_yoke!, rand(8); diffeq)
    # comp_state(y) = (length(y) == 2) ? rand(8) : y; 
    # psys = ProjectedDynamicalSystem(df, [1,2], comp_state)
    yg = range(-10, 10; length = 501)
    grid = ntuple(x -> yg, 8)
    mapper = AttractorsViaRecurrences(df, grid; Δt = 0.1, maximum_iterations = 1e6, Ttr = 400)
    xg = range(-2, 2,length = res)
    yg = range(-2.,2.,length = res)
    grid = (xg, yg)
    bsn = @showprogress [mapper(rand(8)*0.2) for x in xg, y in yg]  
    att = extract_attractors(mapper)
    return @strdict(bsn, att, grid, res)
end


res = 10
params = @strdict res
print_fig(params, "bell_yoke", compute_bell_yoke; ylab= L"\dot{\theta}", xlab= L"\theta", force = true)
att = get_att(params, "bell_yoke", compute_bell_yoke)

diffeq = (reltol = 1e-8,  alg = Vern9(),)
df = CoupledODEs(bell_yoke!, rand(8); diffeq)
y,t = trajectory(df, 500, rand(8)*0.01)
