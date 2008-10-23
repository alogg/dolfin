from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

# Copyright (C) 2005-2007 Anders Logg.
# Licensed under the GNU LGPL Version 2.1.
#
# First added:  2005
# Last changed: 2007-05-14
#
# The bilinear form a(v, U) and linear form L(v) for
# Poisson's equation.
#
# Compile this form with FFC: ffc -l dolfin Poisson.form.

element = FiniteElement("Lagrange", "triangle", 2)

v = TestFunction(element)
u = TrialFunction(element)
f = Function(element)
g = Function(element)

a = dot(grad(v), grad(u))*dx
L = v*f*dx

compile([a, L, M, element], "PoissonP2", {'language': 'dolfin', 'blas': False, 'form_postfix': True, 'precision': '15', 'cpp optimize': False, 'split_implementation': False, 'quadrature_points': False, 'output_dir': '.', 'representation': 'tensor', 'cache_dir': None, 'optimize': False})
