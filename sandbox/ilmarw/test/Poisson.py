from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

element = FiniteElement("Lagrange", "tetrahedron", 1)

v = TestFunction(element)
u = TrialFunction(element)

a = dot(grad(v), grad(u))*dx

compile([a, L, M, element], "Poisson", {'language': 'dolfin', 'blas': False, 'form_postfix': True, 'precision': '15', 'cpp optimize': False, 'split_implementation': False, 'quadrature_points': False, 'output_dir': '.', 'representation': 'tensor', 'cache_dir': None, 'optimize': False})
