using DrWatson
@quickactivate 
using Colors
using ColorSchemes

cmap = ColorScheme([RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )
include("newton_basins.jl")
print_fig(600, 600, cmap, 3, 2000) 
get_Sb(3, 1000)

include("ott_basins.jl")
cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
print_fig(600, 500, cmap, 1200) 
@show get_Sb(1000)

include("duffing_basins.jl")
cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
d = 0.1; F=0.1; ω=0.1;  # smooth boundary
print_fig(600, 500, cmap, d, F, ω, 600) 
@show get_Sb(d, F, ω, 1000)
d = 0.4; F=0.1; ω=0.1;  # smooth boundary
print_fig(600, 500, cmap, d, F, ω, 600) 
@show get_Sb(d, F, ω, 1000)

include("pendulum_basins.jl")
d = 0.2; F = 1.3636363636363635; ω = 0.5 # Parameters for Riddled Basins
cmap = ColorScheme([RGB(0,0,0), RGB(1,1,1)] )
print_fig(600, 500, cmap, d, F, ω, 1000)
@show get_Sb(d, F, ω, 1000)
d = 0.2; F = 1.66; ω = 1. # Parameters for Wada Basins
cmap = ColorScheme([RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)] )
print_fig(600, 500, cmap, d, F, ω, 1000)
@show get_Sb(d, F, ω, 1000)

include("cuencas_hh.jl")
print_fig(600, 600, 0.25, 800) 
print_fig(600, 600, 0.2, 800) 

include("4d_kicked.jl")

include("kicked_rotor.jl")

include("map_feudel.jl")

include("map_grebogi.jl")

include("josephson_junction.jl")

include("magnetic_pendulum.jl")

include("rikitake.jl")

include("thomas.jl")

include("kur_halekotte.jl")


include("lorenz84.jl")
