from dolfin import *

from sys import path 
path.append("../poisson")

from Krylov import *

class Prec: 
    def __init__(self): 
        pass
    def __mul__(self, other): 
        return other.copy()

# Create mesh and finite element
N = 32 
mesh = UnitSquare(N,N)
mesh.disp()
element = FiniteElement("Lagrange", "triangle", 1)
DG = FiniteElement("DG", "triangle", 0)
vector_element = VectorElement("Lagrange", "triangle", 1)
epsilon = 1.0
w_value = 1.0

# Source term
class Source(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)
    def eval(self, values, x):
        values[0] = 0 

# Source term
class BC(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)
    def eval(self, values, x):
        y = x[1]
        values[0] = exp(w_value*y/epsilon) - y*exp(w_value/epsilon) 
        

# Velocity term
class W(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)
    def eval(self, values, x):
        values[0] = 0
        values[1] = w_value  

    def rank(self):
        return 1

    def dim(self, i):
        return 2

# Sub domain for Dirichlet boundary condition
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool(on_boundary)

# Define variational problem
v = TestFunction(element)
w = W(vector_element, mesh)
u = TrialFunction(element)
f = Source(element, mesh)
h = 1.0/N
tau = Function(DG, mesh, 1.0*h**2) 
eps = Function(DG, mesh, epsilon) 

a = dot(w, grad(u))*v*dx + eps*dot(grad(v), grad(u))*dx + tau*dot(dot(w, grad(u)), dot(w, grad(v)))*dx
L = v*f*dx 

# Assemble matrices
A = assemble(a, mesh)
b = assemble(L, mesh)

# Define boundary condition
boundary = DirichletBoundary()
bc_func = BC(element, mesh)
bc = DirichletBC(bc_func, mesh, boundary)
bc.apply(A, b, a)

# create solution vector (also used as start vector) 
x = b.copy()
x.zero()

B = Prec()

# solve the system
regularization_parameter = 1.0/10
#x = Richardson(A, x, b, regularization_parameter, 10e-6, True, 10000)
x = precRichardson(B, A, x, b, regularization_parameter, 10e-6, True, 10000)

# plot the solution
U = Function(element, mesh, x)
plot(U)
interactive()

# Save solution to file
#file = File("conv-diff.pvd")
#file << U


