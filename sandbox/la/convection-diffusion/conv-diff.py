
from dolfin import *

from sys import path 
path.append("../poisson")


from Krylov import *
from BlockLinearAlgebra import *

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
#        values[0] = x[0]
        values[0] = (1- exp((y-1)/epsilon))/(1 - exp(-2/epsilon))
        

# Velocity term
class W(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)
    def eval(self, values, x):
        values[0] = 0
        values[1] = 1  

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

a = dot(w, grad(u))*v*dx + dot(grad(v), grad(u))*dx + tau*dot(dot(w, grad(u)), dot(w, grad(v)))*dx
#a = eps*dot(grad(v), grad(u))*dx 
L = v*f*dx 

# Assemble matrices
A = assemble(a, mesh)
b = assemble(L, mesh)


# Define boundary condition
boundary = DirichletBoundary()
bc_func = BC(element, mesh)
bc = DirichletBC(bc_func, mesh, boundary)
bc.apply(A, b, a)

b.disp()

# create solution vector (also used as start vector) 
x = b.copy()
x.zero()

# create block system
AA = BlockMatrix(1,1); AA[0,0] = A
xx = BlockVector(1);   xx[0]   = x 
bb = BlockVector(1);   bb[0]   = b 

#solve(A,x,b)

BB = Prec()

# solve the system
xx = BiCGStab(AA, xx, bb, 10e-12, True, 1000)
#xx = precondBiCGStab(BB, AA, xx, bb, 10e-12, True, 1000)
#xx = CGN_AA(AA, xx, bb, 10e-12, True)

# plot the solution
U = Function(element, mesh, x)
plot(U)
interactive()

# Save solution to file
#file = File("conv-diff.pvd")
#file << U


