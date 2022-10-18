# BASINS FOR MAJUMDAR-PAPAPETROU DOUBLE BLACK HOLE SYSTEM WITH SYMPLECTIC ALGORITHM*/
function double_black_hole_de!(du, u, p, t)
	
    z1 = p[1]; z2 = p[2]; M1 = p[3]; M2 = p[4]; ∂H∂pϕ = p[5] 
    z, pz, ρ, pρ, ϕ = v
    t_ρ_z1 = (ρ^2+(z-z1)^2)
    t_ρ_z2 = (ρ^2+(z-z2)^2)
    ∂H∂z = pz/(1 + M1/√t_ρ_z1 + M2/√t_ρ_z2)^2
    ∂H∂pz = -2*(1 + M1/√t_ρ_z1 + M2/√t_ρ_z2)*((z-z1)*M1/(t_ρ_z1)^1.5 + (z-z2)*M2/(t_ρ_z2)^1.5)
    ∂H∂ρ = y1/(1 + M1/√t_ρ_z1 + M2/√t_ρ_z2)^2
    ∂H∂pρ = -2*ρ*(1 + M1/√t_ρ_z1 + M2/√t_ρ_z2)*(M1/(t_ρ_z1^1.5) + M2/(t_ρ_z2^1.5)) + pϕ^2/ρ^3/(1 + M1/√t_ρ_z1 + M2/√t_ρ_z2)^2
    ∂H∂ϕ = pϕ/ρ^2/(1 + M1/√t_ρ_z1 + M2/√t_ρ_z2)^2

    dv[1] = ∂H∂pz
    dv[2] = -∂H∂z
    dv[3] = ∂H∂pρ
    dv[4] = -∂H∂ρ
    dv[5] = ∂H∂ϕ
    # dv[6] = ∂H∂pϕ
end

function salida_ρ_z( ρ,  z,  p,  eps) 
	    sal =0 
	if (√(ρ^2+z^2)>10) # escaping photons
		sal=0;
	elseif (abs(z-z1)+abs(ρ)<eps) # upper shadow
		sal=1;
	elseif (abs(z-z2)+abs(ρ)<eps) # lower shadow
		sal=-1;
	else 
		sal=-10;
	end    
return sal
end



# eps=a/10;



# Variables to store trajectories and energy evolution

y0 = range(ρi, ρf, npuntos) 
yp = range(zi, zf, npuntos) 
        
ϕ = 0

U = (1+M1/√t_ρ_z1 + M2/√t_ρ_z2)

if (U^4-pϕ^2/ρ^2<0)
    cuenca_local[l-1][k-1] = -2;
end
p = √(U^4-pϕ^2/ρ^2); # momento total p=√(p_ρ^2+p_z^2)
pρ = p/√(1+(ρ-√3/2)^2/z^2); # tangential shooting from ρ=√(3)/2, z=0 (maximum of the potential) to determine p_ρ
p_z = √(p^2-pρ^2); # p_z from energy conservation

if  z > 0
    pρ = -pρ;
end  
if  ρ - √3/2 < 0
    p_z = -p_z
end
       
     
            
