[![DOI](https://zenodo.org/badge/1024331632.svg)](https://doi.org/10.5281/zenodo.18594252)
## Vector field theory in motion: Revealing latent potentials in football dynamics
This repository contains the MATLAB source code associated with the manuscript "Vector field theory in motion: Revealing latent potentials in soccer dynamics".

This project provides a computational framework to model the collective movement of football players as a vector field, allowing for the reconstruction of a scalar potential field that drives the dynamics of the match.
The code allows the transformation of discrete player trajectories into a continuous field representation to analyze structural properties of team organization.

The framework requires user-provided spatially binned average velocity matrices (`Wx`, `Wy`) as input. Using these data, the main script performs the following steps:
- Visualization of the vector field
- Computation of the rotational
- Verification of Gauss's theorem
- Numerical solution of the scalar potential field
- Visualization of the resulting potential surface

An artificially generated vector field is included in this repository to illustrate the complete workflow. 
All figures displayed in this README are based on this synthetic vector field and are intended for illustrative purposes only.

![vector_field.pdf](https://github.com/user-attachments/files/25185813/vector_field.pdf)


