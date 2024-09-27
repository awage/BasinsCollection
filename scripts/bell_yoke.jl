using DrWatson
@quickactivate
using OrdinaryDiffEq
using Attractors
using CairoMakie
using LaTeXStrings
using ColorSchemes

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
A = 14.; F = 80.; ω = 7. 
Mt(ϕ,dϕ) = abs(ϕ) ≤ π/A ? F*r*sign(dϕ)*cos(ω*ϕ) : 0.

ϕ1, ϕ2, xs, ys, dϕ1, dϕ2, dxs, dys  = u

du[1] = dϕ1
du[2] = dϕ2
du[3] = dxs
du[4] = dys

# (M+m)*ddxs + α*cos(ϕ1)*ddϕ1 + m*l*cos(ϕ2)*ddϕ2 = α*dϕ1^2*sin(ϕ1) + m*l*dϕ2^2*sin(ϕ2) - kx*xs - dx*dxs

# -β*sin(ϕ1)*ddϕ1 - m*l*sin(ϕ2)*ddϕ2 + (M+m)*ddys = β*cos(ϕ1)*dϕ1^2 + m*l*cos(ϕ2)*dϕ2^2 + g*(m + M) - ky*ys - dy*dys

# m*l*cos(ϕ2)*ddx_s + m*l*lc*cos(ϕ1-ϕ2)*ddϕ1 + BC*ddϕ2 - m*l*sin(ϕ2)*ddy_s = m*l*lc*dϕ1^2*sin(ϕ1-ϕ2) - m*g*l*sin(ϕ2) - dc*(dϕ2 - dϕ1)

# α*cos(ϕ1)*ddxs + γ*ddϕ1 + m*l*lc*cos(ϕ1-ϕ2)*ddϕ2 - β*sin(ϕ1)*ddys = -m*l*lc*sin(ϕ1-ϕ2)*dϕ2^2 - α*g*sin(ϕ1) - db*dϕ1 - dc*(dϕ1 - dϕ2) - dbf*2/pi*atan(10^5*dϕ1) - dcf*2/pi*atan(10^5*(dϕ1 - dϕ2)) + Mt

AA = [
α*cos(ϕ1) m*l*cos(ϕ2) (M+m) 0; 
-β*sin(ϕ1) -m*l*sin(ϕ2) 0 M+m; 
m*l*lc*cos(ϕ1-ϕ2) BC m*l*cos(ϕ2) m*l*sin(ϕ2); 
γ m*l*lc*cos(ϕ1-ϕ2) α*cos(ϕ1) -β*sin(ϕ1)]


b1 = α*dϕ1^2*sin(ϕ1) + m*l*dϕ2^2*sin(ϕ2) - kx*xs - dx*dxs; 
b2 = β*cos(ϕ1)*dϕ1^2 + m*l*cos(ϕ2)*dϕ2^2 + g*(m + M) - ky*ys  -  dy*dys;
b3 = m*l*lc*dϕ1^2*sin(ϕ1 - ϕ2) - m*g*l*sin(ϕ2) - dc*(dϕ2 - dϕ1);
b4 = -m*l*lc*sin(ϕ1 - ϕ2)*dϕ2^2 - α*g*sin(ϕ1) - db*dϕ1 - dc*(dϕ1 - dϕ2) - dbf*2/pi*atan(10^5*dϕ1) - dcf*2/pi*atan(10^5*(dϕ1 - dϕ2))+Mt(ϕ1,dϕ1)

du[5:8] = AA\[b1;b2;b3;b4]

end


function compute_bell_yoke(di::Dict)
    @unpack res = di
    diffeq = (reltol = 1e-8,  alg = Vern9(),)
    df = CoupledODEs(bell_yoke!, rand(8); diffeq)
    xg = range(-4.,4.,length = 10001)
    yg = range(-4.,4.,length = 10001)
    psys = ProjectedDynamicalSystem(df, [1,5], zeros(6))
    mapper = AttractorsViaRecurrences(psys, (xg, yg); Δt = 0.1)
    xg = range(-0.1,0.1,length = res)
    yg = range(-1.,1.,length = res)
    grid = (xg, yg)
    bsn, att = basins_of_attraction(mapper, grid)
    return @strdict(bsn, att, grid, res)
end


res = 10
params = @strdict res
print_fig(params, "bell_yoke", compute_bell_yoke; ylab= L"\dot{\theta}", xlab= L"\theta", force = false)
# att = get_att(params, "gear_rattle", compute_gear_rattle)

# diffeq = (reltol = 1e-8,  alg = Vern9(),)
# df = CoupledODEs(bell_yoke!, rand(8); diffeq)
# y,t = trajectory(df, 500, rand(8)*0.01)
