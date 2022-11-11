using DrWatson
@quickactivate 
using Colors
using ColorSchemes

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
