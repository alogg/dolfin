from ufl import *
from ufl.log import set_level
from ffc.compiler.compiler import compile

# Set debug level
set_level(20)

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

# Copyright (C) 2005-2007 Anders Logg.
# Licensed under the GNU LGPL Version 2.1.
#
# First added:  2005
# Last changed: 2008-12-31
#
# The bilinear form a(v, u) and linear form L(v) for
# Poisson's equation.
#
# Compile this form with FFC: ffc -l dolfin Poisson.ufl.

#element = FiniteElement("Lagrange", triangle, 1)
element = VectorElement("Lagrange", triangle, 1)

v = TestFunction(element)
u = TrialFunction(element)
f = Function(element)

a = (inner(grad(v), grad(u)) + inner(v, u))*dx
L = inner(v, f)*dx

compile([a, L, M, element], "PoissonP1", {'log_level': 20, 'format': 'dolfin', 'form_postfix': True, 'precision': '15', 'quadrature_order': 'auto', 'cpp optimize': False, 'cache_dir': None, 'split': False, 'representation': 'auto', 'optimize': False, 'quadrature_rule': None, 'output_dir': '.'}, globals())
