using DrWatson
@quickactivate 
using Colors
using ColorSchemes

include(srcdir("print_fig.jl"))

res = 1200

include("newton_basins.jl")

include("ott_basins.jl")

include("duffing_basins.jl")

include("pendulum_basins.jl")

include("para_pendulum.jl")

include("cuencas_hh.jl")

include("4d_kicked.jl")

include("kicked_rotor.jl")

include("kicked_ott.jl")

include("map_feudel.jl")

include("map_grebogi.jl")

include("josephson_junction.jl")

include("magnetic_pendulum.jl")

include("rikitake.jl")

include("thomas.jl")

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

include("map_cbic.jl")

include("map_kapitaniak.jl")

include("map_cpled_logstc.jl")

include("sprott_memristive.jl")
 
include("coupled_lorenz.jl")

include("map_bogdanov.jl")

include("cyclic_competition.jl")

include("dsnm.jl")

include("split_ring_resonator.jl")

include("eckhardt_9D.jl")

include("kur_halekotte.jl")

include("alfven_complexity.jl")

include("co2_modulated_laser.jl")

include("pumped_laser_dynamics.jl")
 

include("TLNs_network.jl")

include("kuramoto_basins_with_tentacles.jl")

include("map_ashwin.jl")

include("photonic_coupler.jl")

include("megastability_sprott.jl")

include("matryoshka.jl")

include("megastability_sprott.jl")

include("map_mira.jl")

include("map_nash_eq.jl")

include("map_shrimali.jl")

include("map_vcsd_chen.jl")

include("map_bairstow.jl")

include("carpet_oscillator.jl")

include("chimera_reduced_model.jl")

include("gear_rattle.jl")

include("economic_geographic_model.jl")

include("map_lai.jl")
 
include("map_shrimali.jl")

include("decroly_biorhythm.jl")

include("neuron_synapse.jl")

include("souza_impacting.jl")

include("lorenz_computational_chaos.jl")

include("grazing2.jl")

include("bell_yoke2.jl")

include("map_nash_eq.jl")
