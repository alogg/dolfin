from dolfin import *

from sys import path 
path.append("../poisson")


from PyTrilinos import Epetra, EpetraExt, ML, AztecOO 

from Krylov import *

class MLPreconditioner: 
    def __init__(self, prec): 
        self.prec = prec

    def __mul__(self, b):
        x = b.copy()
        x.zero()
        self.prec.Apply(b.vec(),x.vec())

        print "b inner ", b.inner(b) 
        print "x inner ", x.inner(x) 
        return x
    

# Create mesh and finite element
mesh = UnitSquare(2,2)
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



#backend = EpetraFactory.instance()
backend = EpetraFactory.instance()
#backend = PETScFactory.instance()

# Assemble matrices
#A = assemble(a, mesh, None, None, None, None, None, None, None, backend)
A = assemble(a, mesh, backend=backend)

A.disp()


#file=File("A.m") ; file << A;  

#b = assemble(L, mesh, backend=backend) 
#print "b inner ", b.inner(b)

# Define boundary condition
#u0 = Function(mesh, 0.0)
#boundary = DirichletBoundary()
#bc = DirichletBC(u0, mesh, boundary)
#bc.apply(A, b, a)

# create solution vector (also used as start vector) 
#x = b.copy()
#x.zero()

#x = BiCGStab(A, x, b, 10e-12, True, 1000)





# sets up the parameters for ML using a python dictionary
MLList = {
  "max levels"        : 3, 
  "output"            : 10,
  "smoother: type"    : "symmetric Gauss-Seidel",
  "aggregation: type" : "Uncoupled"
}

# creates the preconditioner and computes it
print type(A)
print type(A.mat())
print dir(A.mat())
Prec = ML.MultiLevelPreconditioner(A.mat(), False)
Prec.SetParameterList(MLList)
print "computing ML prec" 
Prec.ComputePreconditioner()
print "done computing ML prec" 
B = MLPreconditioner(Prec)
B.__mul__(x.vec(),b.vec())



dabla = """
# sets up the solver, specifies Prec as preconditioner, and
# solves using CG.
#Solver = AztecOO.AztecOO(A.mat(), x.vec(), b.vec())
#Solver.SetPrecOperator(Prec)
#Solver.SetAztecOption(AztecOO.AZ_solver, AztecOO.AZ_cg)
#Solver.SetAztecOption(AztecOO.AZ_output, 16)
#Solver.Iterate(1550, 1e-5)




regularization_parameter = 1.0/10
x = precRichardson(B, A, x, b, regularization_parameter, 10e-6, True, 10)

r = b-A*x
print "r inner ", r.inner(r) 

#print x 

#A.disp()

# plot the solution
U = Function(element, mesh, x)
plot(U)
interactive()

# Save solution to file
#file = File("poisson.pvd")
#file << U

"""


