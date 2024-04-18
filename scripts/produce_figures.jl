using DrWatson
@quickactivate 
using Colors
using ColorSchemes

include(srcdir("print_fig.jl"))

include("newton_basins.jl")

include("ott_basins.jl")

include("duffing_basins.jl")

include("pendulum_basins.jl")

include("cuencas_hh.jl")

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

include("slim_fractals.jl")

include("sporadical_map.jl")

include("map_henon.jl")

include("cold_atoms.jl")

include("competition_ecology_plancton.jl")

include("open_billiard_sinai.jl")

include("aguirre_billiard.jl")

include("kicked_ott.jl")

include("poon_billiard.jl")

include("hidden_chua.jl")

include("li_sprott.jl")

include("sommerer_basins.jl")

include("black_holes.jl")

include("predator_prey_discrete.jl")

include("lebovitz_mariotti.jl")

include("eckhardt_9D.jl")

include("map_cbic.jl")

include("map_kapitaniak.jl")

include("map_cpled_logstc.jl")

include("sprott_memristive.jl")
 
include("coupled_lorenz.jl")

include("map_bogdanov.jl")

include("cyclic_competition.jl")
