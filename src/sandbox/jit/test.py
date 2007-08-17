from dolfin import *

element = FiniteElement("Lagrange", "triangle", 1)
mesh = UnitSquare(2, 2)

f = Function(element, mesh, 1.0)

v = TestFunction(element)
u = TrialFunction(element)
a = dot(grad(v), grad(u))*dx
L = v*f*dx

A = assemble(a, mesh)
b = assemble(L, mesh)
