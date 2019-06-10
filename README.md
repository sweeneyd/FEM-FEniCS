# FEM-FEniCS
Projects and models I've created in the [FEniCS](https://fenicsproject.org/) finite element framework. Code is from FEniCS 2016.2.0 to 2019.1.0, in which time significant changes have occurred with new dependencies on `pybind11`; earlier code may not run as written in more recent FEniCS builds. Projects include:

- **Poisson**: Documentation for solving the 2D Poisson equation using both a continuous Galerkin finite element method.

- **ParticleTracing**: A mini-project containing a simulation of pressure-driven flow through a pipe and using the solution to trace particles, tied together with Bash scripts.

- **Poission-Nernst-Planck**: A mini-project recreating figures from a 2017 paper by Gao & He and a 2018 paper by Gao & Sun that solves the transient PNP equations for drift-diffusion using a continuous Galerkin and a mixed element formulation (CG/RT), respectively.
