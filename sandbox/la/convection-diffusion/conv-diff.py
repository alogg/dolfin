from PyTrilinos import Epetra, AztecOO, TriUtils, ML 
from dolfin import *

from sys import path 
path.append("../poisson")

from Krylov import *

class MLPreconditioner: 
    def __init__(self, A): 
        # Sets up the parameters for ML using a python dictionary
        MLList = {
              "max levels"        : 5, 
              "output"            : 10,
              "smoother: type"    : "ML symmetric Gauss-Seidel",
              "aggregation: type" : "Uncoupled",
              "ML validate parameter list" : False
        }
        ml_prec = ML.MultiLevelPreconditioner(A.mat(), False)
        ml_prec.SetParameterList(MLList)
        ml_prec.ComputePreconditioner()
        self.ml_prec = ml_prec

    def __mul__(self, b):
        x = b.copy()
        x.zero()
        err = self.ml_prec.ApplyInverse(b.vec(),x.vec())
        if not err == 0: 
            print "err ", err
            return -1 
        return x

# Create mesh and finite element
N = 10 
mesh = UnitSquare(N,N)
element = FiniteElement("Lagrange", "triangle", 1)
DG = FiniteElement("DG", "triangle", 0)
vector_element = VectorElement("Lagrange", "triangle", 1)
epsilon = 1.0/100
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
        c = 1/epsilon
        values[0] = exp(c*y)/exp(c) 

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
tau = Function(DG, mesh, 5.0*h**2) 
eps = Function(DG, mesh, epsilon) 

a = eps*dot(grad(v), grad(u))*dx - dot(w, grad(u))*v*dx + tau*dot(dot(w, grad(u)), dot(w, grad(v)))*dx
L = v*f*dx 


# Assemble matrix and right hand side
backend = EpetraFactory.instance()
A = assemble(a, mesh, backend=backend)
b = assemble(L, mesh, backend=backend)

# Define boundary condition
boundary = DirichletBoundary()
bc_func = BC(element, mesh)
bc = DirichletBC(bc_func, mesh, boundary)
bc.apply(A, b, a)

# create solution vector (also used as start vector) 
x = b.copy()
x.zero()

# create preconditioner
B = MLPreconditioner(A)

# solve the system
regularization_parameter = 1.0/2
x = precRichardson(B, A, x, b, regularization_parameter, 10e-8, True, 20)
#x = precondBiCGStab(B, A, x, b, 10e-6, True, 20)

print x.norm(linf)


# plot the solution
U = Function(element, mesh, x)
plot(U)
interactive()

# Save solution to file
#file = File("conv-diff.pvd")
#file << U


