#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 20:07:54 2019

REFERENCES: https://fenicsproject.org/docs/dolfin/1.6.0/python/demo/documented/mixed-poisson-dual/python/documentation.html
            https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1004.html
            https://fenicsproject.org/pub/tutorial/html/._ftut1008.html
            
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
CG_elem = FiniteElement("CG", mesh.ufl_cell(), 1)
DRT_elem = FiniteElement("RT", mesh.ufl_cell(), 2)

# Assemble mixed element space
W_elem = MixedElement([DRT_elem, DRT_elem, DRT_elem, 
                       CG_elem, CG_elem, CG_elem,
                       R_elem])
W2_elem = MixedElement([DRT_elem, CG_elem, R_elem])
W = FunctionSpace(mesh, W_elem)
W2 = FunctionSpace(mesh, W2_elem)
V = FunctionSpace(mesh, CG_elem)
Q = FunctionSpace(mesh, DRT_elem)

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

bc = DirichletBC(W2.sub(0), Constant((0,0)), boundary)

(sigma0, u0, c0) = TrialFunctions(W2)
(tau_u, v_u, v_c) = TestFunctions(W2)

p0 = interpolate(p_init, V)
n0 = interpolate(n_init, V)

F0 = (dot(tau_u, sigma0) + dot(u0, div(tau_u)))*r*dx
F0 += -(v_u*eps*div(sigma0) + v_u*F*(p0 - n0))*r*dx
F0 += (v_u*c0 + v_c*u0)*r*dx

a0 = lhs(F0)
L0 = rhs(F0)

w0 = Function(W2)
solve(a0 == L0, w0, bc)

(sigma0, u0, c0) = w0.split()

bc = [DirichletBC(W.sub(0), Constant((0,0)), boundary),
      DirichletBC(W.sub(1), Constant((0,0)), boundary),
      DirichletBC(W.sub(2), Constant((0,0)), boundary)]

## Define trial and test functions
(sigma, J_n, J_p, u, n, p, c) = TrialFunctions(W)
(tau_u, tau_n, tau_p, v_u, v_n, v_p, v_c) = TestFunctions(W)

w0 = Function(W)
w = Function(W)

Jp = project(D*grad(p0) + mu*p0*sigma0, Q)
Jn = project(D*grad(n0) - mu*n0*sigma0, Q)

assign(w0.sub(0), sigma0)
assign(w0.sub(1), Jn)
assign(w0.sub(2), Jp)
assign(w0.sub(3), u0)
assign(w0.sub(4), n0)
assign(w0.sub(5), p0)
(sigma0, J_n0, J_p0, u0, n0, p0, c0) = w0.split()

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

a = (inner(J_p, tau_p) + inner(D*p, div(tau_p)))*r*dx 
a += (k*v_p*p - v_p*div(J_p))*r*dx

a += (inner(J_n, tau_n) + inner(D*n, div(tau_n)))*r*dx 
a += (k*v_n*n - v_n*div(J_n))*r*dx

a += (dot(tau_u, sigma) + dot(u, div(tau_u)))*r*dx
a += (-v_u*eps*div(sigma) - v_u*F*(p - n))*r*dx

# Adjustments for Neumann problem
a += (v_u*c + v_c*u)*r*dx

L = f1*v_p*r*dx + f2*v_n*r*dx \
    - inner(mu*n0*sigma0, tau_n)*r*dx \
    + inner(mu*p0*sigma0, tau_p)*r*dx \
    + k*v_p*p0*r*dx \
    + k*v_n*n0*r*dx \
    
keys = ['sigma', 'J_n', 'J_p', 'u', 'n', 'p', 'c']
hfiles = []
for id, solute in enumerate(keys):
    ofile = f'Simulation Results/{solute}/{solute}.pvd'
    hfiles.append(File(ofile))

# Compute solution
t = 0.0
T = 100*dt

storeFiles(hfiles, w0, t)

while t < T:
    solve(a == L, w0, bc)
    t += dt
    storeFiles(hfiles, w0, t)

z = w0.split(True)


mag = lambda vec: project(sqrt(dot(vec,vec)), V)

plt.figure(num=0, figsize=(12,6))
plt.subplot(2,3,1)
p1 = plot(z[5])
cbar = plt.colorbar(p1)
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$P_h$')

plt.subplot(2,3,2)
p2 = plot(z[4])
cbar = plt.colorbar(p2)
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$N_h$')

plt.subplot(2,3,3)
p3 = plot(z[3])
cbar = plt.colorbar(p3)
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$\Psi_h$')

plt.subplot(2,3,4)
#p4 = plot(mag(z[2]), cmap='jet')
#cbar = plt.colorbar(p4)
plot(-z[2])
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$-(\mathbf{J}_p)_h$')

plt.subplot(2,3,5)
#p5 = plot(mag(z[1]), cmap='jet')
#cbar = plt.colorbar(p5)
plot(-z[1])
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$-(\mathbf{J}_n)_h$')

plt.subplot(2,3,6)
#p6 = plot(mag(z[0]), cmap='jet')
#cbar = plt.colorbar(p6)
plot(-z[0])
plt.xlabel('x')
plt.ylabel('y')
plt.title(r'$-\mathbf{\sigma}_h$')

plt.tight_layout()
plt.savefig('Fig3.png', dpi=300)
