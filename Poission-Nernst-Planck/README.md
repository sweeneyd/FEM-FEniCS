# FEM-FEniCS/Poisson-Nernst-Planck

## Project Description
This mini-project contains scripts to to recreate the results from [Gao & He's 2017 paper](https://link.springer.com/article/10.1007/s10915-017-0400-4) and [Gao & Sun's 2018 paper](https://link.springer.com/article/10.1007/s10915-018-0727-5) in which they develop a linearized conservative finite element methods for the Poisson–Nernst–Planck Equations using the [FEniCS 2019.1.0](https://fenicsproject.org/) finite element framework. Specifically, the code mostly recreates Figure 5 from [1] and Figures 3, 4, and 5 from [2] using the using the CG/RT element pairing Anaconda Python 3.7 distribution.

[1] Gao, H., & He, D. (2017). Linearized conservative finite element  methods for the Nernst–Planck–Poisson equations. Journal of Scientific Computing, 72(3), 1269-1289.

[2] Gao, H., & Sun, P. (2018). A Linearized Local Conservative Mixed Finite Element Method for Poisson–Nernst–Planck Equations. Journal of Scientific Computing, 77(2), 793-817.

## File Summary
- `pnp_mwe_CG.py`: Python code with linearized transient PNP solver using a CG method based on [1]
- `pnp_mwe_RT.py`: Python code with linearized transient PNP solver using a CG/RT mixed method based on [2]

## Output Files (`pnp_mwe_RT.py`)
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
- `Fig3_RT.png`: PNG figure showing the solutions for u, n, p, sigma, J_n, and J_p, at time T (300 dpi)

## Output Files (`pnp_mwe_CG.py`)
- `c/c.pvd`: VTK-formatting file for solution for undetermined constant for Neumann problem (c)
- `c/cXXXXXX.vtu`: VTK-formatting file for solution for undetermined constant for Neumann problem (c)
- `u/u.pvd`: VTK-formatting file for solution for scalar potential field (u)
- `u/uXXXXXX.vtu`: VTK-formatting file for solution for scalar potential field (u)
- `n/n.pvd`: VTK-formatting file for solution for negatively-charged species (n)
- `n/nXXXXXX.vtu`: VTK-formatting file for solution for negatively-charged species (n)
- `p/p.pvd`: VTK-formatting file for solution for positively-charged species (p)
- `p/pXXXXXX.vtu`: VTK-formatting file for solution for positively-charged species (p)
- `Fig3_CG.png`: PNG figure showing the solutions for u, n, and p at time T (300 dpi)
