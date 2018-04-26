# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:10:19 2017

@author: Dan

To run: 1.) Install FEniCS in Docker (via Anaconda forge)
        2.) Export .msh file from GMSH with NONE of the export options checked
        3.) In the directory containing the .msh file 
                >> source activate fenicsproject
        4.) dolfin-convert square_hole2.msh square_hole2.xml

Great FEniCS Resource: https://fenicsproject.org/pub/course/lectures/2015-08-logg-beijing/
"""

from dolfin import *
import matplotlib.pyplot as plt
import numpy as np

## Import Mesh ---------------------------------------------------------------
ifilename = 'square_hole2.msh'
ofilename = ifilename.split('.')[0]+'.xml'
mesh = Mesh(ofilename)

## Get boundaries/subdomain
# This pulls out the labels of the subdomains and boundaries from the auxiliary files
subdomains = MeshFunction("size_t", mesh, "%s_physical_region.xml"%ofilename.split('.')[0])
boundaries = MeshFunction("size_t", mesh, "%s_facet_region.xml"%ofilename.split('.')[0])

## Initialize solution space and parameters-----------------------------------
# Define function space (i.e. 'CG' := Continuous Galerkin)
V = FunctionSpace(mesh, 'CG', 1)

# Define test/trial functions functions
u = TrialFunction(V)
v = TestFunction(V)

# Define integration labels
dx = Measure('dx', domain=mesh, subdomain_data=subdomains)
ds = Measure('ds', domain=mesh, subdomain_data=boundaries)

# Define constants
D = Constant(1)
alpha = Constant(1e-3)
g = Constant(1e-4)
cL = Constant(20e-6)
cR = Constant(0)

## Define boundary conditions ------------------------------------------------
boundary_conditions = {5: {'Neumann': g},
                       2: {'Dirichlet': cR},
                       4: {'Dirichlet': cL}}
DBCs = []
integrals_N = []
integrals_R_a = []
integrals_R_L = []
for i in boundary_conditions:
    if 'Dirichlet' in boundary_conditions[i]:
        bc = DirichletBC(V, boundary_conditions[i]['Dirichlet'], boundaries, i)
        DBCs.append(bc)

    if 'Neumann' in boundary_conditions[i]:
        if boundary_conditions[i]['Neumann'] != 0:
            g = boundary_conditions[i]['Neumann']
            integrals_N.append(g*v*ds(i))

## Define the weak form  -----------------------------------------------------
a = dot(grad(u), grad(v))*dx - alpha*u*v*dx + sum(integrals_R_a)
L = sum(integrals_N) + sum(integrals_R_L)

## Solve system  -------------------------------------------------------------
# Define a new function to store the solution
u = Function(V)

# Solve the system and store into u, applying the Dirichlet B.C.'s
solve(a==L, u, DBCs)

## Save solution to file in VTK format
vtkfile = File('Solution/example2D_steadyState.pvd')
vtkfile << u

## Make pretty plots  --------------------------------------------------------
# Plot mesh
fig = plt.figure(num=2)
plot(mesh)
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$p=1$ trianglar mesh')
plt.savefig('plotMesh.png', dip=300)
plt.savefig('plotMesh.eps')

# Plot solution
fig = plt.figure(num=1)
p = plot(u, cmap='jet')
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$D\nabla^2u - \alpha u = 0$')
plt.colorbar(p, label=r'$u(x,y)$')
plt.savefig('plotSolution.png', dip=300)
plt.savefig('plotSolution.eps')

## Unsteady state problem  ---------------------------------------------------
# Initialize time-stepping using Backward Euler's Method
u = TrialFunction(V)
v = TestFunction(V)
u0 = interpolate(Constant(0), V)
time = np.linspace(0, 1, 11)
theta = 1/2
a = u*v*dx + dt*(D*dot(grad(u), grad(v))*dx - alpha*u*v*dx) 
A = assemble(a)
L = u0*v*dx + dt*g*v*ds(5)

u1 = Function(V)
vtkfile = File('Solution/example2D_CrankNicholson.pvd')
T = 1
dt = 1e-2
t = 0
while t <= T:
    b = assemble(L)
    for i in DBCs:
        i.apply(A, b)
    solve(A, u1.vector(), b)
    t += dt
    vtkfile << (u1, t)
    u0.assign(u1)   