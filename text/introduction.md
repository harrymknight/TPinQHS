# Modelling Thermopower in Quantum Hall Systems

_Harry Knight_

The quantum Hall effect, discovered in 1980 by von Klitzing et al., manifests itself as emergent plateaus of Hall conductance as a function of magnetic field or electron density in two-dimensional electron systems (2DES) [(1)](https://en.wikipedia.org/wiki/Quantum_Hall_effect). Systems which display this effect are termed quantum Hall systems (QHS) and are of interest due to the complex electronic distributions underlying them. The simplest of these systems, and the focus of this dashboard, are integer quantum Hall systems (IQHS).

An IQHS is generated when a 2DES, either an inversion layer of a MOSFET [(2)](https://en.wikipedia.org/wiki/Depletion_region#Depletion_width_in_MOS_capacitor) or at the interface between layers of a heterojunction [(3)](https://en.wikipedia.org/wiki/Heterojunction), is kept at extremely low temperatures (less than 100mK) and exposed to a strong magnetic field (typically greater than 10T) perpendicular to its surface (see directly below). The resulting electronic distribution is determined by single particle quantum states whose wavefunctions, satisfy the Schrödinger equation [(4)](https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation), and therefore the 2nd order ODE,

$$
E \psi_n = \frac{\hbar \omega_c}{2} \bigg{\lbrack}-{l_B}^2 {\partial_x}^2 + \frac{{(k{l_B}^2 + x)}^2}{{l_B}^2} \bigg{\rbrack} \psi_n
$$

with $n$ denoting the energy level of the wavefunction and $\hbar \omega_c$ and $l_B$ respectively the characteristic energy and length of the system.

The solutions of the 2nd order ODE, called eigenfunctions [(5)](https://en.wikipedia.org/wiki/Eigenfunction), are the eigenfunctions of a quantum harmonic oscillator [(6)](https://en.wikipedia.org/wiki/Quantum_harmonic_oscillator) displaced by the cycloctron orbit centre of the wavefunction $x_k = -k {l_B}^2$ with $k$ the wave number [(7)](https://en.wikipedia.org/wiki/Wavenumber) of the wavefunction. 

The resulting energy levels, called Landau levels [(8)](https://en.wikipedia.org/wiki/Landau_quantization) are therefore said to be degenerate [(9)](https://en.wikipedia.org/wiki/Degenerate_energy_levels), meaning an energy level is common to multiple quantum states $\Leftrightarrow$ wavefunctions. As an edge of the system is approached the boundary conditions [(10)](https://en.wikipedia.org/wiki/Boundary_value_problem) of the 2nd order ODE change leading to quadratic growth of the Landau levels. This is presented in the energy spectrum ($E_{n{,}k}$ vs $x_k$ graph) seen below which is a cross-section of the system.