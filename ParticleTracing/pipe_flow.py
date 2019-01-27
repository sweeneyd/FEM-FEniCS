from dolfin import *
import sys

# Load mesh from file
ifilename = sys.argv[1]
ofilename = ifilename.split('.')[0]+'.xml'
mesh = Mesh(ofilename)

## Get boundaries/subdomain
boundaries = MeshFunction("size_t", mesh, "%s_facet_region.xml"%ofilename.split('.')[0])

# Define function spaces
V = VectorElement("CG", mesh.ufl_cell(), 2)
Q = FiniteElement("CG", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, V*Q)

# No-slip boundary condition for velocity
noslip = Constant((0, 0, 0))
bc0 = DirichletBC(W.sub(0), noslip, boundaries, 1)

# Inflow boundary condition for velocity
# x0 = 1
inflow = Constant(("0", "0.0", "-1.0"))
bc1 = DirichletBC(W.sub(0), inflow, boundaries, 3)

# Boundary condition for pressure at outflow
# x0 = 0
zero = Constant(0)
bc2 = DirichletBC(W.sub(1), zero, boundaries, 2)

# Collect boundary conditions
bcs = [bc0, bc1, bc2]

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
f = Constant((0, 0, 0))
a = (inner(grad(u), grad(v)) - div(v)*p + q*div(u))*dx
L = inner(f, v)*dx

# Compute solution
w = Function(W)
solve(a == L, w, bcs)

# # Split the mixed solution using a shallow copy
(u, p) = w.split()

save_dir = "Results/"

# Save solution in .h5 format
Uoutput_file = HDF5File(mesh.mpi_comm(), "velocity.h5", "w")
Uoutput_file.write(u, "velocity")
Uoutput_file.close()

Poutput_file = HDF5File(mesh.mpi_comm(), "presure.h5", "w")
Poutput_file.write(p, "pressure")
Poutput_file.close()
