from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

# Copyright (c) 2005-2007 Anders Logg
# Licensed under the GNU LGPL Version 2.1
#
# First added:  2005
# Last changed: 2007-04-30
#
# The bilinear form a(v, u) and Linear form L(v) for the Stokes
# equations using a mixed formulation (equal-order stabilized).
#
# Compile this form with FFC: ffc -l dolfin Stokes.form

vector = VectorElement("Lagrange", "triangle", 1)
scalar = FiniteElement("Lagrange", "triangle", 1)
system = vector + scalar

(v, q) = TestFunctions(system)
(u, p) = TrialFunctions(system)

f = Function(vector)
h = Function(scalar)

beta  = 0.2
delta = beta*h*h

a = (dot(grad(v), grad(u)) - div(v)*p + q*div(u) + delta*dot(grad(q), grad(p)))*dx
L = dot(v + mult(delta, grad(q)), f)*dx

compile([a, L, M, element], "Stokes", options={'language': 'dolfin', 'blas': False, 'form_postfix': True, 'precision': '15', 'cpp optimize': False, 'split_implementation': False, 'quadrature_points': False, 'output_dir': '.', 'representation': 'tensor', 'cache_dir': None, 'optimize': False}, global_variables=globals())
