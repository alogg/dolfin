"""This demo program solves the mixed formulation of Poisson's
equation:

    sigma + grad(u) = 0
         div(sigma) = f

The corresponding weak (variational problem)

    <tau, sigma> - <div(tau), u> = 0       for all tau
                 <w, div(sigma)> = <w, f>  for all w

is solved using BDM (Brezzi-Douglas-Marini) elements of degree q (tau,
sigma) and DG (discontinuous Galerkin) elements of degree q - 1 for
(w, u).

Original implementation: ../cpp/main.cpp by Anders Logg and Marie Rognes
"""

__author__    = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__      = "2007-11-14 -- 2008-12-19"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__   = "GNU LGPL Version 2.1"

from dolfin import *

# Create mesh and define function spaces
mesh = UnitSquare(16, 16)

# Works!
BDM = FunctionSpace(mesh, "BDM", 1)
DG = FunctionSpace(mesh, "DG", 0)

# Does not work!
#BDM = FunctionSpace(mesh, "BDM", 1)
#DG = FunctionSpace(mesh, "CG", 1)

# Does not work!
#BDM = VectorFunctionSpace(mesh, "CG", 2)
#DG = FunctionSpace(mesh, "DG", 0)

# Does not work!
BDM = VectorFunctionSpace(mesh, "CG", 2)
DG = FunctionSpace(mesh, "CG", 1)

V = BDM + DG

# Define variational problem
(tau, w) = TestFunctions(V)
(sigma, u) = TrialFunctions(V)
f = Expression("500.0*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", V = DG)

#a = (dot(tau, sigma) - div(tau)*u + w*div(sigma))*dx
#L = w*f*dx

a = dot(tau, sigma)*dx + dot(tau, grad(u))*dx - dot(grad(w), sigma)*dx + 0.00001*w*u*dx
L = w*f*dx

#a = (dot(tau, sigma) + inner(tau, grad(u)) - inner(grad(w), sigma))*dx #+ 0.01*w*u*dx
#L = w*f*dx

#zero = Constant(mesh, 0)
#bc = DirichletBC(V.sub(1), zero, DomainBoundary())

# Compute solution
problem = VariationalProblem(a, L)
(sigma, u) = problem.solve().split()

# Project sigma for post-processing
P1 = VectorFunctionSpace(mesh, "CG", 1)
sigma_proj = project(sigma, P1)
# FIXME: Strange result here when interpolating!

print "||div(sigma)|| =", norm(div(sigma))
print "||f|| =", norm(f)
print "||div(sigma) - f|| =", norm(div(sigma) - f)

V = VectorFunctionSpace(mesh, "CG", 2)

print norm(div(project(sigma, V)) - f)

# Plot solution
plot(sigma_proj)
interactive()
plot(u)
interactive()

# Save solution to pvd format
f3 = File("sigma.pvd")
f4 = File("u.pvd")
f3 << sigma_proj
f4 << u
