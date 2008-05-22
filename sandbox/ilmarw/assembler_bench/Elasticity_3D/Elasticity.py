from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

# Copyright (c) 2005 Johan Jansson (johanjan@math.chalmers.se)
# Licensed under the GNU LGPL Version 2.1
#
# Modified by Anders Logg 2006-2007
#
# First added:  2005
# Last changed: 2007-04-18
#
# The bilinear form for classical linear elasticity (Navier).
# Compile this form with FFC: ffc -l dolfin Elasticity.form.

element = VectorElement("Lagrange", "tetrahedron", 1)

v = TestFunction(element)
u = TrialFunction(element)
f = Function(element)

E  = 10.0
nu = 0.3

mu    = E / (2*(1 + nu))
lmbda = E*nu / ((1 + nu)*(1 - 2*nu))

def epsilon(v):
    return 0.5*(grad(v) + transp(grad(v)))

def sigma(v):
    return 2*mu*epsilon(v) + lmbda*mult(trace(epsilon(v)), Identity(len(v)))

a = dot(grad(v), sigma(u))*dx
L = dot(v, f)*dx

compile([a, L, M, element], "Elasticity", "tensor", "dolfin", {'quadrature_points=': False, 'blas': False, 'precision=': '15', 'optimize': False})
