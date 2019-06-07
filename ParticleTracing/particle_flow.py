from dolfin import *
import sys
import numpy as np
from scipy.optimize import fsolve

# Load mesh from file
ifilename = sys.argv[1]
velocity_filename = sys.argv[2]
tracers = int(sys.argv[3])
S = float(sys.argv[4])
delta_s = float(sys.argv[5])

ofilename = ifilename.split('.')[0]+'.xml'
mesh = Mesh(ofilename)

V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)
Q2 = FunctionSpace(mesh, 'CG', 2)
T = TensorFunctionSpace(mesh, "CG", 1)

velocity = Function(V)
velocity_file = HDF5File(mesh.mpi_comm(), velocity_filename, "r")
velocity_file.read(velocity, "velocity")
velocity_file.close()

abs_velocity = project(sqrt(inner(velocity, velocity)), Q2)
normal_velocity = project(velocity/abs_velocity, V)
J = project(grad(velocity), T)

def vel(x):
    point = Point(x)
    return velocity(point)

def jac(x):
    point = Point(x)
    return J(point).reshape((3,3))

def func(x, x_0, dt):
    u = vel(x)
    return x - x_0 - u*dt

def funcjac(x, x_0, dt):
    return np.eye(3)-jac(x)*dt

N = int(S/delta_s)
eps = 0.01
x_min = -1.0
x_max = 1.0
y_min = -1.0
y_max = 1.0
z_min = 0.0
z_max = 10.0


fout = open("particles.csv", 'w')
for tr in range(tracers):
    print("", tr)
    theta = np.random.uniform(0, 2*np.pi)
    radius = np.random.uniform(0, 1.0)
    x_0 = np.array([radius*np.cos(theta),
                    radius*np.sin(theta),
                    eps,])
    x_n = x_0
    t_n = 0.0
    for i in range(N):
        u_n = vel(x_n)
        print(tr, t_n, x_n[0], x_n[1], x_n[2], u_n[0], u_n[1], u_n[2], file=fout, sep=',')
        if x_n[0] >= x_max:
            break
        dt = delta_s/np.linalg.norm(u_n)
        t_n += dt
        # Crank--Nicolson scheme
        x_n1 = x_n + u_n*dt
        x_n2 = fsolve(func, x_n1, args=(x_n, dt), fprime=funcjac)
        x_n = x_n2

        if mesh.bounding_box_tree().compute_first_entity_collision(Point(x_n)) < mesh.num_cells():
        # if x_n[2] <= z_min+eps or x_n[2] >= z_max-eps \
        #    or np.sqrt(x_n[0]**2 + x_n[1]**2) >= np.sqrt(x_max**2 + y_max**2)-eps:
            break
    print("\n\n", file=fout)
fout.close()
