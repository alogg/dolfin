
from PyTrilinos import Epetra, AztecOO, TriUtils, ML 
from dolfin import *

from sys import path 
path.append("../poisson")
from Krylov import *


class MLPreconditioner: 
    def __init__(self, A): 
        # Sets up the parameters for ML using a python dictionary
        MLList = {
              "max levels"        : 3, 
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
mesh = UnitSquare(20,20)
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
#backend = PETScFactory.instance()

# Assemble matrices
A = assemble(a, mesh, backend=backend)
b = assemble(L, mesh, backend=backend) 

#file = File("A.m"); file <<A

# Define boundary condition
u0 = Function(mesh, 0.0)
boundary = DirichletBoundary()
bc = DirichletBC(u0, mesh, boundary)
bc.apply(A, b, a)

# create solution vector (also used as start vector) 
x = b.copy()
x.zero()

B = MLPreconditioner(A)
x = precondconjgrad(B, A, x, b, 10e-6, True, 100)
#x = conjgrad(A, x, b, 10e-6, True, 100)

# plot the solution
U = Function(element, mesh, x)
plot(U)
interactive()

# Save solution to file
file = File("poisson.pvd")
file << U




