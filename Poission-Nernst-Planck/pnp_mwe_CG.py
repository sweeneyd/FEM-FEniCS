#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 11:34:31 2019

** Gao, H., & He, D. (2017). Linearized conservative finite element 
   methods for the Nernst–Planck–Poisson equations. Journal of Scientific 
   Computing, 72(3), 1269-1289.

@author: Dan
"""

from dolfin import *
import matplotlib.pyplot as plt

def storeFiles(hfiles, w_, t):
    _w_temp = w_.split(True)
    for file_no, c_ in enumerate(_w_temp):
        hfiles[file_no] << (c_, t)

## Create mesh
mesh = UnitSquareMesh(16, 16)

def boundary(x, on_boundary):
        return on_boundary

# Define function spaces and mixed (product) space
R_elem = FiniteElement("R", mesh.ufl_cell(), 0)
CG_elem = FiniteElement("CG", mesh.ufl_cell(), 2)

# Assemble mixed element space
W_elem = MixedElement([CG_elem, CG_elem, CG_elem, R_elem])
W2_elem = MixedElement([CG_elem, R_elem])
W = FunctionSpace(mesh, W_elem)
W2 = FunctionSpace(mesh, W2_elem)
V = FunctionSpace(mesh, CG_elem)

# Create material property
p_init = Expression("(x[0] > xcut)*(x[1] > ycut) + err", 
                    xcut=0.75, ycut=11./20., 
                    err=1e-6, degree=CG_elem.degree())
n_init = Expression("(x[0] > xcut)*(x[1] < ycut) + err", 
                    xcut=0.75, ycut=9./20., 
                    err=1e-6, degree=CG_elem.degree())

# Define geometric coefficients
r = Expression("1.0", degree=CG_elem.degree())

## Plot exact solution
#u_exact = Expression("exp(t)*cos(pi*x[0])*cos(pi*x[1]) - t*t*t*cos(2*pi*x[0])", t=1.0, degree=CG_elem.degree())
#p_exact = Expression("2*pi*pi*exp(t)*cos(pi*x[0])*cos(pi*x[1])", t=1.0, degree=CG_elem.degree())
#n_exact = Expression("4*pi*pi*t*t*t*cos(2*pi*x[0])", t=1.0, degree=CG_elem.degree())
#
#U = project(u_exact, V)
#P = project(p_exact, V)
#N = project(n_exact, V)
#
#plt.subplot(1,3,1)
#plot(U)
#
#plt.subplot(1,3,2)
#plot(P)
#
#plt.subplot(1,3,3)
#plot(N)
#
#plt.tight_layout()
#

D = Constant(1.0)
mu = Constant(1.0)
eps = Constant(1.0)
F = Constant(1.0)

D = Constant(1.0)
mu = Constant(1.0)
eps = Constant(1.0)
F = Constant(1.0)

(u0, c0) = TrialFunctions(W2)
(v_u, v_c) = TestFunctions(W2)

p0 = interpolate(p_init, V)
n0 = interpolate(n_init, V)

F0 = dot(grad(v_u), grad(u0))*r*dX \
    - v_u*(p0 - n0)*r*dx \
    + (v_u*c0 + v_c*u0)*r*dx
    

a0 = lhs(F0)
L0 = rhs(F0)

w0 = Function(W2)
solve(a0 == L0, w0, [])

(u0, c0) = w0.split()

## Define trial and test functions
(u, n, p, c) = TrialFunctions(W)
(v_u, v_n, v_p, v_c) = TestFunctions(W)

w0 = Function(W)
w = Function(W)

assign(w0.sub(0), u0)
assign(w0.sub(1), n0)
assign(w0.sub(2), p0)
assign(w0.sub(3), c0)
(u0, n0, p0, c0) = w0.split()

## Solve FEM solution
#f1 = p_exact - div(grad(p_exact) - p*grad(u_exact))
#f2 = n_exact - div(grad(n_exact) + n*grad(u_exact))
f1 = Constant(0.)
f2 = Constant(0.)
#f1 = Expression("sin(10*(1-x[0]*x[1]))", degree=CG_elem.degree())
#f2 = Expression("0.5*sin(10*x[0])*cos(10*x[1])", degree=CG_elem.degree())

dt = 2.0e-3
k = Constant(1.0/dt)

# Define variational form
Fi = (k*v_p*(p - p0) + dot(D*grad(p), grad(v_p)) + dot(mu*p0*grad(u), grad(v_p)))*r*dx

Fi += (k*v_n*(n - n0) + dot(D*grad(n), grad(v_n)) - dot(mu*n0*grad(u), grad(v_n)))*r*dx

Fi += dot(grad(v_u), grad(u))*r*dx \
    - v_u*F/eps*(p - n)*r*dx

# Adjustments for Neumann problem
Fi += (v_u*c + v_c*u)*r*dx

a = lhs(Fi)
L = rhs(Fi)

keys = ['u', 'n', 'p', 'c']
hfiles = []
for id, solute in enumerate(keys):
    ofile = f'Simulation Results/CG/{solute}/{solute}.pvd'
    hfiles.append(File(ofile))

# Compute solution
t = 0.0
T = 100*dt

storeFiles(hfiles, w0, t)

while t < T:
    solve(a == L, w0, [])
    t += dt
    storeFiles(hfiles, w0, t)

z = w0.split(True)


mag = lambda vec: project(sqrt(dot(vec,vec)), V)

plt.figure(num=0, figsize=(12,6))
plt.subplot(2,3,1)
p1 = plot(z[2])
cbar = plt.colorbar(p1)
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$P_h$')

plt.subplot(2,3,2)
p2 = plot(z[1])
cbar = plt.colorbar(p2)
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$N_h$')

plt.subplot(2,3,3)
p3 = plot(z[0])
cbar = plt.colorbar(p3)
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$\Psi_h$')

plt.subplot(2,3,4)
#p4 = plot(mag(z[2]), cmap='jet')
#cbar = plt.colorbar(p4)
plot(-grad(z[2]))
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$-(\mathbf{J}_p)_h$')

plt.subplot(2,3,5)
#p5 = plot(mag(z[1]), cmap='jet')
#cbar = plt.colorbar(p5)
plot(-grad(z[1]))
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$-(\mathbf{J}_n)_h$')

plt.subplot(2,3,6)
#p6 = plot(mag(z[0]), cmap='jet')
#cbar = plt.colorbar(p6)
plot(-grad(z[0]))
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$-\mathbf{\sigma}_h$')

plt.tight_layout()
plt.savefig('Fig3_CG.png', dpi=300)