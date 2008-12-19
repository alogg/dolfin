"""This demo program solves Poisson's equation

    - div grad u(x) = f(x)

on the unit interval with source f given by

    f(x) = 9.0*DOLFIN_PI*DOLFIN_PI*sin(3.0*DOLFIN_PI*x[0]);

and boundary conditions given by

    u(x) = 0 for x = 0
    du/dx = 0 for x = 1
"""

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-11-28 -- 2008-12-19"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU LGPL Version 2.1"

from dolfin import *

# Create mesh and finite element
mesh = UnitInterval(50)

V = FunctionSpace(mesh, "Lagrange", 1)

# Source term
class Source(Function):
    def eval(self, values, x):
        values[0] = 9.0*DOLFIN_PI*DOLFIN_PI*sin(3.0*DOLFIN_PI*x[0])

# Flux term
class Flux(Function):
    def eval(self, values, x):
        values[0] = 3.0*DOLFIN_PI*cos(3.0*DOLFIN_PI*x[0])

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool(on_boundary) and bool(x[0] < DOLFIN_EPS)

# Define variational problem
v = TestFunction(V)
u = TrialFunction(V)
f = Source(V)
g = Flux(V)

a = dot(grad(v), grad(u))*dx
L = v*f*dx + v*g*ds

# Define boundary condition
u0 = Constant(mesh, 0.0)
bc = DirichletBC(V, u0, DirichletBoundary())

# Solve PDE and plot solution
pde = LinearPDE(a, L, bc)
u = pde.solve()

# Save solution to file
file = File("poisson.pvd")
file << u

# Plot solution
plot(u, interactive=True, rescale=False)
