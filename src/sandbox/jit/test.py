from dolfin import *

element = FiniteElement("Lagrange", "triangle", 1)

v = TestFunction(element)
u = TrialFunction(element)
a = dot(grad(v), grad(u))*dx

mesh = UnitSquare(2, 2)
A = assemble(a, mesh)
