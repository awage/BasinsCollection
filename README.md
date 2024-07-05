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
* Riddled basins from: https://doi.org/10.1016/0167-2789(94)90047-7
* Duffing oscillator basins 
* Forced pendulum basins
* Henon Heiles open Hamiltonian escape basins 
* 4D kicked map rotor (1998) `4d_kicked.jl`: https://doi.org/10.1016/S0960-0779(97)00058-1
* 4D kicked map rotor from: https://doi.org/10.1103/PhysRevE.54.71
* 4D kicked map rotor from: https://doi.org/10.1016/0167-2789(87)90108-4
* A map with Wada: https://doi.org/10.1103/PhysRevE.58.3060
* Map with fractal boundary from: https://doi.org/10.1016/0375-9601(83)90945-3
* Fractal boundary in the josephson junction: https://doi.org/10.1103/PhysRevA.36.2455 
* Fractal Boundary from the magnetic pendulum: https://doi.org/10.1103/PhysRevLett.111.194101
* Rikitake oscillator basins, unpublished.
* Thomas cyclical oscillator basins, unpublished
* Kuramoto oscillators on the UK power grid: https://doi.org/10.1088/2632-072X/ac080f
* Lorenz 84 system: https://doi.org/10.1063/1.2953589
* Slim fractals: https://doi.org/10.1103/PhysRevX.7.021040
* Sporadical map: https://doi.org/10.1103/PhysRevLett.82.3597
* Hénon map: https://doi.org/10.1016/j.physrep.2014.02.007
* Cold atoms: https://doi.org/10.1103/PhysRevA.95.013629
* Competition ecology plancton: https://doi.org/10.1086/319929
* Open Sinai billiard: https://doi.org/10.1103/PhysRevA.38.930
* Aguirre billiard: https://doi.org/10.1103/PhysRevE.67.056201
* Poon billiard: https://doi.org/10.1142/S0218127496000035
* Hidden attractors in Chua oscillators: https://doi.org/10.1016/j.physleta.2011.04.037
* Li Sprott oscillator: https://doi.org/10.1142/S0218127416502333 
* Riddled Basins from: Sommerer, John C. "The end of classical determinism." Johns Hopkins APL Technical Digest 16.4 (1995): 333.
* Binary Black Holes escape basins: https://doi.org/10.1103/PhysRevD.98.084050
* Discrete predator Prey system: https://doi.org/10.1016/j.chaos.2022.112833
* Extreme multistable system: https://doi.org/10.1103/PhysRevE.85.035202
* Disipative nontwist map: https://doi.org/10.48550/arXiv.2211.06921
* 6D Shear Flow model: https://doi.org/10.1017/jfm.2013.38
* 9D model of Fluid dynamics: https://doi.org/10.1088/1367-2630/6/1/056 and https://doi.org/10.1103/PhysRevE.91.052903 (this last paper has a basins that I cannot reproduce)
* Wada basins in the cubic 2D map (`map_cbic.jl`): http://dx.doi.org/10.1016/j.physleta.2013.03.027
* Riddled basins of discrete 2D system: https://doi.org/10.1103/PhysRevE.57.R6253
* Coupled logistic maps with riddled basins (Fig 15 of the paper has an error): https://doi.org/10.1103/PhysRevE.57.2713
* Sprott-memristor model: https://doi.org/10.1016/j.chaos.2022.111834
* Parametrically forced pendulum: http://dx.doi.org/10.1142/S0218127411030167
* Bogdanov map:  https://doi.org/10.1142/S021812749300074X
* Coupled Lorenz systems:  https://doi.org/10.1103/PhysRevE.96.062203
* Rock-Paper-Scisors competition model `cyclic_competition.j`:  https://doi.org/10.1063/1.5045366
* Split Ring Resonator model `split_ring_resonator.jl`:  https://doi.org/10.1063/5.0157489
* Alfven Complexity `alfven_complexity.jl`: # ALFVÉN BOUNDARY CRISIS: https://doi.org/10.1142/S0218127402005303

