from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

element = FiniteElement("Lagrange", "tetrahedron", 1)

v = TestFunction(element)
u = TrialFunction(element)

a = dot(grad(v), grad(u))*dx

compile([a, L, M, element], "Poisson", "tensor", "dolfin", {'quadrature_points=': False, 'blas': False, 'precision=': '15', 'optimize': False})
