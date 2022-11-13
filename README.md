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
* Riddled basins from E. Ott, et al. I Physica D 76 (1994) 384-410
* Duffing oscillator basins 
* Forced pendulum basins
* Henon Heiles open Hamiltonian escape basins 
* 4D kicked map rotor: Feudel, U., Grebogi, C., Poon, L., & Yorke, J. A. (1998). Dynamical properties of a simple mechanical system with a large number of coexisting periodic attractors. Chaos, Solitons & Fractals, 9(1-2), 171-180.
* 4D kicked rotor from: Map with more than 100 coexisting low-period periodic attractors, Ulrike Feudel,  Celso Grebogi, Brian R. Hunt, and James A. Yorke
* A map with Wada:  Basin bifurcation in quasiperiodically forced systems Ulrike Feudel, Annette Witt, Ying-Cheng Lai, and Celso Grebogi PRE 28, 1998
* Map with fractal boundary from: C. Grebogi, S. W. McDonald, E. Ott, J. A. Yorke, Final state sensitivity: An obstruction to predictability, Physics Letters A, 99, 9, 1983
* Fractal boundary in the josephson junction: Cell-mapping techniques applied to the rf-driven Josephson junction, Physical review A, Atomic, molecular, and optical physics · October 1987  DOI: 10.1103/PhysRevA.36.2455 
* Fractal Boundary from the magnetic pendulum
* Rikitake oscillator basins, unpublished.
* Thomas cyclical oscillator basins, unpublished
* Kuramoto oscillators on the UK power grid:  Halekotte, Feudel. 
* Lorenz 84 system:  Multistability, phase diagrams, and intransitivity in the Lorenz-84 low-order atmospheric circulation model Chaos 18, 033121 (2008); https://doi.org/10.1063/1.2953589
* Slim fractals: Phys. Rev. X 7, 021040 (2017) Slim Fractals: The Geometry of Doubly Transient Chaos Xiaowen Chen,1, ∗ Takashi Nishikawa,1, 2 and Adilson E. Motter
* Sporadical map: Alexandre, Ott et. al. 
* Hénon map 
* Cold atoms
* Competition ecology plancton: vol. 157, no. 5 the american naturalist, Fundamental Unpredictability in Multispecies Competition, Jef Huisman, and Franz J. Weissing
* Open Sinai billiard: PHYSICAL REVIEW A VOLUME 38, NUMBER 2 JULY 15, 1988 Fractal boundaries for exit in Hamiltonian Dynamics, Siegfried Bleher, Celso Grebogi, Edward Ott, and Reggie Brown
* Aguirre billiard: PHYSICAL REVIEW E 67, 056201, Limit of small exits in open Hamiltonian systems. Jacobo Aguirre and Miguel A. F. Sanjuán
* Kicked rotor:  Multi dimensioned intertwined basin boundaries: basin structure of THE KICKED DOUBLE ROTOR  Celso GREBOGI a, Eric KOSTELICH b, Edward OTT a,c and James A. YORKE b  Physica 25D 1987
* Poon billiard: Poon, Campos, Ott, Grebogi IJBC, vol6. nº2 1996 251-265
* Hidden attractors in Chua oscillators:  Physics Letters A Volume 375, Issue 23, 6 June 2011, Pages 2230-2233 Localization of hidden Chuaʼs attractors A.Leonov N.V. Kuznetsov  V.I.Vagaitsev
* Li Sprott oscillator:  International Journal of Bifurcation and Chaos, Vol. 26, No. 14 (2016) 1650233 (11 pages) DOI: 10.1142/S0218127416502333 Crisis in Amplitude Control Hides in Multistability Int. J. Bifurcation Chaos 2016.26.  Chunbiao Li, Julien Clinton Sprott, Hongyan Xing. 

