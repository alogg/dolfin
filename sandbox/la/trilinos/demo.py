from dolfin import *

from sys import path 
path.append("../poisson")



from Krylov import *

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

# Neumann boundary condition
class Flux(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)
    def eval(self, values, x):
        if x[0] > DOLFIN_EPS:
            values[0] = 25.0*sin(5.0*DOLFIN_PI*x[1])
        else:
            values[0] = 0.0

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool(on_boundary and x[0] < DOLFIN_EPS)

# Define variational problem
v = TestFunction(element)
u = TrialFunction(element)
f = Source(element, mesh)
g = Flux(element, mesh)

a = dot(grad(v), grad(u))*dx
L = v*f*dx + v*g*ds



backend = EpetraFactory.instance()

# Assemble matrices
A = assemble(a, mesh, backend)
b = assemble(L, mesh, backend)


# Define boundary condition
u0 = Function(mesh, 0.0)
boundary = DirichletBoundary()
bc = DirichletBC(u0, mesh, boundary)
bc.apply(A, b, a)

# create solution vector (also used as start vector) 
x = b.copy()
x.zero()

# create block system
#AA = BlockMatrix(1,1); AA[0,0] = A
#xx = BlockVector(1);   xx[0]   = x 
#bb = BlockVector(1);   bb[0]   = b 


# solve the system
#xx = BiCGStab(AA, xx, bb, 10e-12, True, 1000)


x = BiCGStab(A, x, b, 10e-12, True, 1000)

print x 

# plot the solution
U = Function(element, mesh, x)
plot(U)
interactive()

# Save solution to file
file = File("poisson.pvd")
file << U


