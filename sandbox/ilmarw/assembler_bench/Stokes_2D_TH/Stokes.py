from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

# Copyright (c) 2005-2007 Anders Logg
# Licensed under the GNU LGPL Version 2.1
#
# First added:  2005
# Last changed: 2007-04-23
#
# The bilinear form a(v, u) and Linear form L(v) for the Stokes
# equations using a mixed formulation (Taylor-Hood elements).
#
# Compile this form with FFC: ffc -l dolfin Stokes.form

P2 = VectorElement("Lagrange", "triangle", 2)
P1 = FiniteElement("Lagrange", "triangle", 1)
TH = P2 + P1

(v, q) = TestFunctions(TH)
(u, p) = TrialFunctions(TH)

f = Function(P2)

a = (dot(grad(v), grad(u)) - div(v)*p + q*div(u))*dx
L = dot(v, f)*dx

compile([a, L, M, element], "Stokes", "tensor", "dolfin", {'quadrature_points=': False, 'blas': False, 'precision=': '15', 'optimize': False})
