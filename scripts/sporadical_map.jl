using DynamicalSystems

function Ms(x)
    y=0.
     if x ≤ 0; y=9*x/(4-5*x); end
     if 0 < x ≤ 4/9; y=9/4*x; end
     if 4/9 < x ≤ 5/9; y=81/4*(x-x^2)-4; end
     if  x > 5/9; y=9/4*(1-x); end
     return y
 end

# Basin bifurcation in quasiperiodically forced systems Ulrike Feudel, Annette Witt, Ying-Cheng Lai, and Celso Grebogi PRE 28, 1998
function chaotic_map(dz,z, p, n)
    xn = z[1]; yn = z[2]
    λ=p[1];
    dz[1]=Ms(xn)
    dz[2]=λ*yn+sin(2*π*xn)
    return
end

# dummy function to keep the initializator happy
function chaotic_map_J(J,z0, p, n)
   return
end
ds = DiscreteDynamicalSystem(chaotic_map,[1., 0.], [1.1] , chaotic_map_J)
integ  = integrator(ds)

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


function print_fig(w, h, cmap, res)
    xg=range(0.,1.,length=600)
    yg=range(-2.,6.,length=600)
    integ.p[1] = 1.1
    bsn = [get_color(integ,x,y) for x in xg, y in yg]

    fig = Figure(resolution = (w, h))
    ax = Axis(fig[1,1], ylabel = L"$\dot{x}$", xlabel = L"x", yticklabelsize = 30, 
            xticklabelsize = 30, 
            ylabelsize = 30, 
            xlabelsize = 30, 
            xticklabelfont = "cmr10", 
            yticklabelfont = "cmr10")
    if isnothing(cmap)
        heatmap!(ax, xg, yg, bsn, rasterize = 1)
        else
            heatmap!(ax, xg, yg, bsn, rasterize = 1, colormap = cmap)
        end
    save(string(projectdir(), "/plots/sporadical_",res,".png"),fig)
end



print_fig(600,600, nothing, 500) 
