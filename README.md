# BasinsCollection

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> BasinsCollection

It is authored by Wagemakers.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add(DrWatson) # install globally, for using `quickactivate`
   julia> Pkg.activate(path/to/this/project)
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.


# Content: 

* Basins for the Newton root finding algorithm

   ![](./plots/newton__res=1200.png)

* Riddled basins `ott_basins.jl`: https://doi.org/10.1016/0167-2789(94)90047-7

   ![basins_riddle_ott](./plots/basins_riddle_ott_res=1200.png)

* Riddled Basins from `sommerer_basins.jl`: Sommerer, John C. "The end of classical determinism." Johns Hopkins APL Technical Digest 16.4 (1995): 333.
   
   ![basin_sommerer](./plots/basin_sommerer_res=1200.png)

* Duffing oscillator basins `duffing_basins.jl`

   ![duffing_0.05](./plots/duffing_0.05_F=0.098_d=0.05_res=1200_ω=1.15.png)

* Forced pendulum basins `pendulum_basins.jl`

   ![pendulum](./plots/pendulum_F=1.66_d=0.2_res=1200_ω=1.0.png)

* Hénon map `map_henon.jl`: https://doi.org/10.1016/j.physrep.2014.02.007

   ![henon](./plots/henon_j=0.9_res=1200_μ=1.08.png)

* Henon Heiles open Hamiltonian escape basins `cuencas_hh.jl`

   ![hh](./plots/hh_E=0.25_res=1200.png)

* 4D kicked map rotor (1998) `4d_kicked.jl`: https://doi.org/10.1016/S0960-0779(97)00058-1

   ![kicked_rotor_4d](./plots/kicked_rotor_4d_res=1200_ρ=6.5.png)

* A map with hundred of attractors: `kicked_rotor.jl`: https://doi.org/10.1103/PhysRevE.54.71

   ![kicked_rotor](./plots/kicked_rotor_f0=4.0_res=1200_ν=0.02.png)

* 4D kicked map rotor from `kicked_ott.jl`: https://doi.org/10.1016/0167-2789(87)90108-4

   ![kicked_rotor_ott](./plots/kicked_rotor_ott_f0=0.1_res=1200_ν=0.1.png)

* A map with Wada `map_feudel.jl`: https://doi.org/10.1103/PhysRevE.58.3060

   ![feudel](./plots/feudel_a=0.0015_r=3.83_res=1200.png)

* Map with fractal boundary `map_grebogi.jl`: https://doi.org/10.1016/0375-9601(83)90945-3

   ![grebogi](./plots/grebogi_res=1200.png)

* Fractal boundary in the josephson junction `josephson_junction.jl`: https://doi.org/10.1103/PhysRevA.36.2455 

   ![josephson](./plots/josephson_idc=1.88_irf=10.2_res=1200_Ω=1_β=25.png)

* Fractal Boundary from the magnetic pendulum `magnetic_pendulum.jl`: https://doi.org/10.1103/PhysRevLett.111.194101

   ![mag_pend](./plots/mag_pend_N=3_d=0.3_res=1200_α=0.2_γ=1_ω=0.5.png)

* Rikitake oscillator basins `rikitake.jl`, unpublished.

   ![rikitake](./plots/rikitake_res=1200_α=1.0_μ=0.5.png)

* Thomas cyclical oscillator basins `thomas.jl`, unpublished

   ![thomas](./plots/thomas_b=0.166_res=1200.png)

* Kuramoto oscillators on the UK power grid `kur_halekotte.jl`: https://doi.org/10.1088/2632-072X/ac080f

   ![basins_kur](./plots/basins_kur_ni=14_res=300.png)

* Lorenz 84 system `lorenz84.jl`: https://doi.org/10.1063/1.2953589

   ![lorenz84](./plots/lorenz84_F=6.85_G=1.29_a=0.25_b=4.0_res=1200.png)

* Slim fractals `slim_fractals.jl`: https://doi.org/10.1103/PhysRevX.7.021040

   ![slim_fractal](./plots/slim_fractal_res=1200_μ=0.2.png)

* Sporadical map `sporadical_map.jl`: https://doi.org/10.1103/PhysRevLett.82.3597

   ![sporadical](./plots/sporadical_res=1200.png)

* Cold atoms `cold_atoms.jl`: https://doi.org/10.1103/PhysRevA.95.013629

   ![cold_atoms](./plots/cold_atoms_res=1200_vx=0.1_x=-500.png)

* Competition ecology plancton `competition_ecology_plancton.jl`: https://doi.org/10.1086/319929

   ![comptetion_model](./plots/comptetion_model_res=1200.png)

* Open Sinai billiard `open_billiard_sinai.jl`: https://doi.org/10.1103/PhysRevA.38.930

   ![open_sinai](./plots/open_sinai_res=1200_Δ=0.8.png)

* Aguirre billiard `aguirre_billiard.jl`: https://doi.org/10.1103/PhysRevE.67.056201

   ![open_disks](./plots/open_disks_res=1200_ww=0.001.png)

* Poon billiard `poon_billiard.jl`: https://doi.org/10.1142/S0218127496000035

   ![poon_disks](./plots/poon_disks_res=1200_ww=0.1.png)

* Hidden attractors in Chua oscillators `hidden_chua.jl`: https://doi.org/10.1016/j.physleta.2011.04.037

   ![hidden_chua](./plots/hidden_chua_m0=-0.177_m1=-1.15_res=1200_α=8.46_β=12.1_γ=0.00516.png)

* Li Sprott oscillator `li_sprott.jl`: https://doi.org/10.1142/S0218127416502333 

   ![li_sprott](./plots/li_sprott_a=13.0_b=0.55_res=1200.png)

* Binary Black Holes escape basins `black_holes.jl`: https://doi.org/10.1103/PhysRevD.98.084050

   ![black_holes](./plots/black_holes_res=1200_Δpϕ=0.03.png)

* Discrete predator Prey system `predator_prey_discrete.jl`: https://doi.org/10.1016/j.chaos.2022.112833

   ![pred_prey](./plots/pred_prey_k=7.94_m=0.104_res=1200.png)

* Disipative nontwist map `dsnm.jl`: https://doi.org/10.1103/PhysRevE.107.024216

   ![dsnm](./plots/dsnm_a=0.55_b=0.45_res=2000.png)

* 6D Shear Flow model `lebovitz_mariotti.jl`: https://doi.org/10.1017/jfm.2013.38

   ![lebovitz](./plots/lebovitz_Re=307.0_res=1200.png)

* 9D model of Fluid dynamics `eckhardt_9D.jl`: https://doi.org/10.1088/1367-2630/6/1/056 and https://doi.org/10.1103/PhysRevE.91.052903 (this last paper has a basins that I cannot reproduce)

   ![eckhardt](./plots/eckhardt_Re=425.0_res=1200.png)

* Wada basins in the cubic 2D map (`map_cbic.jl`): http://dx.doi.org/10.1016/j.physleta.2013.03.027
  
  ![cbic_map](./plots/cbic_map_j=0.66_res=1200_μ=2.9.png)

* Riddled basins of discrete 2D system `map_kapitaniak.jl`: https://doi.org/10.1103/PhysRevE.57.R6253

   ![kapitaniak](./plots/kapitaniak_d1=-0.935_d2=-0.935_l=1.41_pp=-1.41_res=1200.png)

* Coupled logistic maps with riddled basins (Fig 15 of the paper has an error) `map_cpled_logstc.jl`:  https://doi.org/10.1103/PhysRevE.57.2713

   ![cplog](./plots/cplog_a=3.6_res=1200_ε=-1.0.png)

* Sprott-memristor model `sprott_memristive.jl`: https://doi.org/10.1016/j.chaos.2022.111834

   ![hidden_memristor](./plots/4d_sprott_memristor_g=0.03_m=11_r=5.8_res=1200_α=1_β=0.05_γ=0.5.png)

* Parametrically forced pendulum `para_pendulum.jl`: http://dx.doi.org/10.1142/S0218127411030167

   ![parametric_pendulum](./plots/parametric_pendulum_A_1=0.0_A_2=4.1_res=1200_ω=1.5.png)

* Bogdanov map `map_bogdanov.jl`:  https://doi.org/10.1142/S021812749300074X

   ![bogdanov](./plots/bogdanov_k=1.2_res=2000_ε=0.0125_μ=-0.1.png)

* Coupled Lorenz systems `coupled_lorenz.jl`:  https://doi.org/10.1103/PhysRevE.96.062203

   ![coupled_lorenz](./plots/coupled_lorenz_res=1200_α=10_β=24.8_γ=2.67_ε=1.1.png)

* Rock-Paper-Scisors competition model `cyclic_competition.j`:  https://doi.org/10.1063/1.5045366

   ![cyc_comp](./plots/cyc_comp_r=3.35_res=1200.png)

* Split Ring Resonator model `split_ring_resonator.jl`:  https://doi.org/10.1063/5.0157489

   ![split_ring_res](./plots/split_ring_res_res=1200_μ=35_ω=1.03.png)

* Alfven Complexity `alfven_complexity.jl`: https://doi.org/10.1142/S0218127402005303

   ![alfven](./plots/alfven_res=1200_ν=0.0175.png)

* Pump modulation erbium doped fiber laser  (basins are slightly different) `pumped_laser_dynamics.jl`: https://doi.org/10.1016/j.physleta.2009.10.061

   ![basins_lasers](./plots/basins_lasers_f=70200.0_m=0.8_res=1200.png)

* CO2 modulated laser model (the basins in the paper are wrong due to numerical instabilities in the model) `co2_modulated_laser.jl`: https://doi.org/10.1063/5.0093727

   ![basins_CO2_laser](./plots/basins_CO2_laser_f=0.005_m=0.02_res=1200.png)

* Threshold-Linear Networks with multistable patterns `TLNs_network.jl`: https://doi.org/10.1371/journal.pone.0264456 

   ![basins_tlns](./plots/basins_tlns_res=1200_θ=1.0.png)

* Basins with tentacles `kuramoto_basins_with_tentacles.jl`:  https://doi.org/10.1103/PhysRevLett.127.194101

   ![basins_tentacles](./plots/basins_tentacles_N=40_res=1200.png)

* Riddled basins in coupled cuadratic map `map_ashwin.jl`: https://doi.org/10.1088/0951-7715/9/3/006

   ![cpashwin](./plots/cpashwin_res=1200_α=0.7_ε=0.5_ν=1.28.png)

* Active Photonic Couplers `photonic_coupler.jl`: https://doi.org/10.1103/PhysRevA.100.043834

   ![photonic_coupler](./plots/photonic_coupler_k=5_res=1200_α=2_β=1.5_ε=0.5.png)

* Multistable Chimera in reduced Kuramoto model (two coupled populations) `chimera_reduced_model.jl`: https://doi.org/10.1088/1367-2630/18/2/022002

   ![chimera_states](./plots/chimera_states_res=1200_β=0.025_μ=0.55_ν=0.45.png)

* Carpet oscillator `carpet_oscillator.jl`: https://doi.org/10.1007/s12043-018-1581-6

   ![carpet](./plots/carpet_res=1200.png)

* Megastability: nested attractors: `megastability_sprott.jl`: https://doi.org/10.1140/epjst/e2017-70037-1

   ![megastable](./plots/megastable_res=1200.png)

* Matryoshka multistability: `matryoshka.jl`: https://doi.org/10.1016/j.chaos.2024.115412

   ![matryoshka](./plots/matryoshka_res=1200.png)

* Bairstow application: `map_bairstow.jl`: https://doi.org/10.1063/1.166158

   ![bairstow_map](./plots/bairstow_map_a=0.8_res=1200.png)

* Quadratic map for basin bifurcations: `map_mira.jl`: https://doi.org/10.1142/S0218127494000241

   ![mira_map](./plots/mira_map_a=-0.42_b=-1.32_res=1200.png)

* Geographic Economic competition model: `economic_geographic_model.jl`: http://dx.doi.org/10.1016/j.matcom.2014.01.004

   ![basins_economic](./plots/basins_economic_res=1500_Φ=0.085.png)

* Multistability in a dynamic Cournot game: `map_nash_eq.jl`: 

   ![nash_equilibrium](./plots/nash_equilibrium_res=1200_λ=0.5_μ=1.95.png)

* Nonsmooth models of gear rattle oscilator: `gear_rattle.jl`: https://doi.org/10.1142/S021812740902283X

   ![gear_rattle](./plots/gear_rattle_res=1200_β=0.6_δ=0.6_ε=0.1.png)

* Piecewise Smooth Dynamical System: `map_lai.jl`: https://doi.org/10.1063/1.2985853

   ![lai_pcw](./plots/lai_pcw_res=1200.png)

* Basin bifurcations in quasiperiodically forced coupled systems: `map_shrimali.jl`: https://doi.org/10.1103/PhysRevE.72.036215

   ![shrimali](./plots/shrimali_res=1200_α=3.25_β=0.01_ϵp=0.5.png)

* Basins of attraction in a gear ratlle oscillator: `gear_rattle_souza.jl`: https://doi.org/10.1177/107754630100700605 

   ![gear_rattle2](./plots/gear_rattle2_res=1200.png)

* Basins of attraction of different ringing schemes of the church bell: `bell_yoke2.jl`: https://doi.org/10.1016/j.ijimpeng.2015.06.008

   ![bell_yoke](./plots/bell_yoke_Tmax=150_lr=-0.03_res=600.png)

* Basins of a multiply regulated biochemical system: `decroly_biorhythm.jl`: https://doi.org/10.1073/pnas.79.22.6917   

   ![biorhythm](./plots/biorhythm_ks=1.99_res=1200.png)

* Multistability in discrete chaotic systems using numerical integration with variable symmetry: `map_vcsd_chen.jl`: https://doi.org/10.1016/j.chaos.2022.112794

  ![chen](./plots/chen_S=0.77_res=1200.png)

* Adaptive synapse-based neuron model with heterogeneous multistability and riddled basins: `neuron_synapse.jl`: https://doi.org/10.1063/5.0125611

   ![neuron_synapse](./plots/neuron_synapse_B=2_c=1.8_g=1.7_res=1200.png)

* Grazing chaos in impacting system: `souza_impacting.jl`: https://doi.org/10.1016/j.chaos.2007.01.022

   ![impacting](./plots/impacting_res=1200.png)

* Simplified discretized Lorenz model: `lorenz_computational_chaos.jl`:  https://doi.org/10.1016/0167-2789(89)90072-9

   ![lorenz_complexity](./plots/lorenz_complexity_res=1200.png)

* Cournot economic game in 2D: `map_bischi_cournot.jl`: https://doi.org/10.1016/S0960-0779(98)00130-1

![cournot2d](./plots/cournot2d_res=1200_μ1=3.55_μ2=3.55.png)

* Chaotic Gyrostat basins: `gyrostat.jl`: https://doi.org/10.3390/math10111914

![gyrostat](./plots/gyrostat=1200.png)
