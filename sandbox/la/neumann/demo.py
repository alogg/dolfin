"""This demo program solves Poisson's equation

    - div grad u(x, y) = f(x, y)

on the unit square with source f given by

    f(x, y) = 500*exp(-((x-0.5)^2 + (y-0.5)^2)/0.02)

and boundary conditions given by

    du/dn(x, y) = 0               
"""

from dolfin import *

# Create mesh and finite element
mesh = UnitSquare(32, 32)
element = FiniteElement("Lagrange", "triangle", 1)

# Source term
class Source(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)
    def eval(self, values, x):
        dx = x[0] - 0.5
        dy = x[1] - 0.5
        values[0] = 500.0*exp(-(dx*dx + dy*dy)/0.02)

# Define variational problem
v = TestFunction(element)
u = TrialFunction(element)
f = Source(element, mesh)

a = dot(grad(v), grad(u))*dx
m = v*u*dx
L = v*f*dx 

x = Vector() 
b = assemble(L, mesh) 
A = assemble(a, mesh) 
M = assemble(m, mesh) 

z = b.copy() 
z.assign(1.0) 
y = b.copy() 
y.assign(0.0) 

#A.mult(b,z)
c  = z.inner(b)/b.size()
y.assign(c)
print "type(b) ", type(b)
print "type(c) ", type(c)
b -= y

A.mult(y, z)

print "zzzzzzzzzzzzzzz"
print z.norm()




solve(A, x, b, gmres, amg)

b.disp()
x.disp()

plot(u)
u  = Function(element, mesh, x)


# Save solution to file
file = File("poisson.pvd")
file << u

# Hold plot
plot(u)
interactive()
