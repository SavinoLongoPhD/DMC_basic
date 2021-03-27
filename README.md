# DMC_basic
A simple Diffusion Monte Carlo program for teaching molecular structure.

It allows to calculated energy and wavefunctions for atoms and molecules with one or two valence electrons.
Examples are:
One-nucleus: H, He, Li+
Two-nuclei: H_2^+, HeH+, H2, LiH
Three-nuclei: H_3^+, LiD_2^+

Released as a teaching tool dfor the class "Spectroscopy and Computer Modeling of Molecular Systems", M.S. Physics, Bari University, Italy.

The program: DMC_Basic.f
Language: fortran77

Details of the algorithm are found in the reference paper Longo, GM, et al. "The unbiased Diffusion Monte Carlo: a versatile tool for two-electron systems confined in different geometries.", EPJ D 2021. Please make reference to the paper if this code is used as a research tool.

The paper is open-access and can be downloaded here https://link.springer.com/article/10.1140/epjd/s10053-021-00095-7

DMC_Basic.f allows to perform many experiments, provided the two-electron system is in its singlet ground state.

A list of quantum systems to study, with the related coulomb energy functions, is provided in the code.
It includes H, H2+, H2.
A simple empty-core pseudo-potential for Li is proposed to make experiments with molecules like LiH+, LiH, Li2+.

For the one-electron case, excited states are readily selected using symmetry-based nodal planes. For details: Longo GM et al, Plasma Sources Science and Technology 24.6 (2015): 065019.

This program uses atomic units

---------------------------------------------------------------------------------------------------------------------

How to use the program: 

To use the program, you must be able to write the explicit expression of the various Coulombian attractions and repulsions in the system you choose to study.

However, a fairly extensive set of examples is found in the program itself.

So you can experiment by simply uncommenting a few lines of your choice.

The variable ndim at the beginning of the program must be set equal to 3 times the number of electrons in the system therefore = 3 for a one-electron system, and = 6 for two electrons. In the case of molecules, a cycle is used, which allows scanning of different internuclear distances "d".
The program prints, during execution, the internuclear distance parameter (which you can possibly ignore if you are not using it) and the estimate of the energy of the atomic-molecular system in atomic units. It is possible, by activating some lines of the program, to create histograms that allow for example to represent the wave function in the system.

Note that this is a Monte Carlo program, consequently the results, even if in principle they represent the exact solution of the wave equation, are affected by statistical errors. You can see them by observing, for example, the trend of the results as a function of the internuclear distance.
To reduce this fluctuation effect, you can use a longer computation time, or average several results. Refer to the publication reported above to understand the calculation method.

