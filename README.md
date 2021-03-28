# DMC_basic
A simple Diffusion Monte Carlo program for teaching molecular structure.

It allows to calculated energy and wavefunctions for atoms and molecules with one or two valence electrons.
Examples are:
One-nucleus: H, He, Li+.
Two-nuclei: H_2^+, HeH+, H2, LiH.
Three-nuclei: H_3^+, LiD_2^+

This version was developed as a teaching tool for the class "Spectroscopy and Computer Modeling of Molecular Systems", M.S. Physics, Bari University, Italy.

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
