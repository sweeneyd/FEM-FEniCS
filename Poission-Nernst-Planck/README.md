# FEM-FEniCS/Poisson-Nernst-Planck

## Project Description
This mini-project contains scripts to to recreate the results from [Gao & Sun's 2018 paper](https://link.springer.com/article/10.1007/s10915-018-0727-5)* in which they develop a linearized local conservative mixed finite element method for the Poisson–Nernst–Planck Equations using the [FEniCS 2019.1.0](https://fenicsproject.org/) finite element framework. Specifically, the code mostly recreates Figures 3, 4, and 5 from this article using the using the CG/RT element pairing Anaconda Python 3.7 distribution.

*Gao, H., & Sun, P. (2018). A Linearized Local Conservative Mixed Finite Element Method for Poisson–Nernst–Planck Equations. Journal of Scientific Computing, 77(2), 793-817.

## File Summary
- `pnp_mwe.py`: Python code with linearized PNP solver

## Output Files
- `c/c.pvd`: VTK-formatting file for solution for undetermined constant for Neumann problem (c)
- `c/cXXXXXX.vtu`: VTK-formatting file for solution for undetermined constant for Neumann problem (c)
- `u/u.pvd`: VTK-formatting file for solution for scalar potential field (u)
- `u/uXXXXXX.vtu`: VTK-formatting file for solution for scalar potential field (u)
- `n/n.pvd`: VTK-formatting file for solution for negatively-charged species (n)
- `n/nXXXXXX.vtu`: VTK-formatting file for solution for negatively-charged species (n)
- `p/p.pvd`: VTK-formatting file for solution for positively-charged species (p)
- `p/pXXXXXX.vtu`: VTK-formatting file for solution for positively-charged species (p)
- `sigma/sigma.pvd`: VTK-formatting file for solution for electrical flux vector (sigma)
- `sigma/sigmaXXXXXX.vtu`: VTK-formatting file for solution for electrical flux vector (sigma)
- `J_n/J_n.pvd`: VTK-formatting file for solution for mass flux vector for species n (J_n)
- `J_n/J_nXXXXXX.vtu`: VTK-formatting file for solution for mass flux vector for species n (J_n)
- `J_p/J_p.pvd`: VTK-formatting file for solution for mass flux vector for species p (J_p)
- `J_p/J_pXXXXXX.vtu`: VTK-formatting file for solution for mass flux vector for species p (J_p)
- `Fig3.png`: PNG figure showing the solutions for u, n, p, sigma, J_n, and J_p, at time T (300 dpi)
