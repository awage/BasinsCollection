using DrWatson
@quickactivate 
using Colors
using ColorSchemes

include(srcdir("print_fig.jl"))

res = 1200

file_arr = [
"newton_basins.jl",

"ott_basins.jl",

"duffing_basins.jl",

"pendulum_basins.jl",

"para_pendulum.jl",

"cuencas_hh.jl",

"4d_kicked.jl",

"kicked_rotor.jl",

"map_feudel.jl",

"map_grebogi.jl",

"josephson_junction.jl",

"magnetic_pendulum.jl",

"rikitake.jl",

"thomas.jl",

"lorenz84.jl",

"slim_fractals.jl",

"sporadical_map.jl",

"map_henon.jl",

"cold_atoms.jl",

"competition_ecology_plancton.jl",

"open_billiard_sinai.jl",

"aguirre_billiard.jl",

"kicked_ott.jl",

"poon_billiard.jl",

"hidden_chua.jl",

"li_sprott.jl",

"sommerer_basins.jl",

"black_holes.jl",

"predator_prey_discrete.jl",

"lebovitz_mariotti.jl",

"map_cbic.jl",

"map_kapitaniak.jl",

"map_cpled_logstc.jl",

"sprott_memristive.jl",
 
"coupled_lorenz.jl",

"map_bogdanov.jl",

"cyclic_competition.jl",

"dsnm.jl",

"eckhardt_9D.jl",

"kur_halekotte.jl"]

Threads.@threads for k in 1:length(file_arr) 
    include(file_arr[k])
end
