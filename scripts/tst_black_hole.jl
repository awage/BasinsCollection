# BASINS FOR MAJUMDAR-PAPAPETROU DOUBLE BLACK HOLE SYSTEM WITH SYMPLECTIC ALGORITHM*/
function double_black_hole_de!(du, u, p, t)
	
    z1 = p[1]; z2 = p[2]; M1 = p[3]; M2 = p[4]; ∂H∂pϕ = p[5] 
    z, pz, ρ, pρ, ϕ = v

    z = pz/(1+M1/sqrt(ρ^2+(z-z1)^2) + M2/sqrt(ρ^2+(z-z2)^2))/(1+M1/sqrt(ρ^2+(z-z1)^2)+M2/sqrt(ρ^2+(z-z2)^2));
    ∂H∂pz = -2*(1+M1/sqrt(ρ^2+(z-z1)^2)+M2/sqrt(ρ^2+(z-z2)^2)) * ((z-z1)*M1/((ρ^2+(z-z1)^2))^1.5+(z-z2)*M2/((ρ^2+(z-z2)^2))^1.5)
    ∂H∂ρ = y1/(1+M1/sqrt(ρ^2+(z-z1)^2) + M2/sqrt(ρ^2+(z-z2)^2)^2)  
    ∂H∂pρ = -2*ρ*(1+M1/sqrt(ρ^2+(z-z1)^2)+M2/sqrt(ρ^2+(z-z2)^2))*(M1/((ρ^2+(z-z1)^2)^1.5)+ M2/((ρ^2+(z-z2)^2)^1.5)) + pϕ^2/ρ^3/(1+M1/sqrt(ρ^2+(z-z1)^2) + M2/sqrt(ρ^2+(z-z2)^2))/(1+M1/sqrt(ρ^2+(z-z1)^2)+M2/sqrt(ρ^2+(z-z2)^2))
    ∂H∂ϕ = pϕ/ρ^2/(1+M1/sqrt(ρ^2+(z-z1)^2)+M2/sqrt(ρ^2+(z-z2)^2))/(1+M1/sqrt(ρ^2+(z-z1)^2)+M2/sqrt(ρ^2+(z-z2)^2));

    dv[1] = ∂H∂pz
    dv[2] = -∂H∂z
    dv[3] = ∂H∂pρ
    dv[4] = -∂H∂ρ
    dv[5] = ∂H∂ϕ
    # dv[6] = ∂H∂pϕ
end

function salida_ρ_z( ρ,  z,  p,  eps) 
	    sal =0 
	if (sqrt(ρ^2+z^2)>10) # escaping photons
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



eps=a/10;

nt=int(tf/h);


# Variables to store trajectories and energy evolution

       y0 = range(ρi, ρf, npuntos) 
       yp = range(zi, zf, npuntos) 
        
        for (k=1;k<=npuntos;k=k+1)
            
        {	cuenca_local[l-1][k-1]=-10;

            ρ=y0[k-1];
	    z=yplocal[l-1];

	    ϕ=0;
	p_ϕ=pϕ;# to avoid contamination from the integration
	U=(1+M1/sqrt(ρ*ρ+(z-z1)*(z-z1))+M2/sqrt(ρ*ρ+(z-z2)*(z-z2)));
        if (pow(U,4)-p_ϕ*p_ϕ/ρ/ρ<0)
        cuenca_local[l-1][k-1]=-2;
        
        p=sqrt(pow(U,4)-p_ϕ*p_ϕ/ρ/ρ); # momento total p=sqrt(p_ρ^2+p_z^2)
        p_ρ=p/sqrt(1+(ρ-sqrt(3)/2)/z*(ρ-sqrt(3)/2)/z); # tangential shooting from ρ=sqrt(3)/2, z=0 (maximum of the potential) to determine p_ρ
        p_z=sqrt(p*p-p_ρ*p_ρ); # p_z from energy conservation
        
        if (z>0)
            p_ρ=-p_ρ;   
        
        if (ρ-sqrt(3)/2<0)
            p_z=-p_z;
       
     
            
