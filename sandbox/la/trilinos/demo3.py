from PyTrilinos import Epetra, AztecOO, TriUtils, ML 
from sys import path 
path.append("../poisson")



from dolfin import *
from BlockLinearAlgebra import *
from Krylov import *
from MinRes import *

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
 
class SaddlePrec: 
    def __init__(self, M, N):  
        self.n = 2 
        self.M = M
        self.N = N
        self.M_prec = MLPreconditioner(M)
        self.N_prec = MLPreconditioner(N)

    def __mul__(self, b):
        if (not isinstance(b, BlockVector)):
            raise TypeError, "Can not multiply DiagBlockMatrix with %s" % (str(type(b)))
        if (not b.n == self.n):
            raise ValueError, "dim(BlockVector) != dim(Prec) "  

        b0 = b.data[0]
        b1 = b.data[1]

        x0 = self.M_prec*b0
        x1 = self.N_prec*b1

        xx = BlockVector(x0, x1)
        return xx 

# Boundary value
class Source(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)

    def eval(self, values, x):
        values[0] = 1.0
        values[1] = 0.0

    def rank(self):
        return 1

    def dim(self, i):
        return 2


N = 1 
mesh = UnitSquare(N,N) 
shape = "triangle"
V = VectorElement("CG", shape, 2)
Q = FiniteElement("CG", shape, 1)

# Test and trial functions
v = TestFunction(V)
u = TrialFunction(V)
q = TestFunction(Q)
p = TrialFunction(Q)
tau = Function(Q, mesh, 0.0)
f = Source(V, mesh)


# velocity terms 
a00 = dot(v, u)*dx + dot(grad(v), grad(u))*dx

# Divergence constraint 
a10 = div(u)*q*dx

# gradient of p 
a01 = div(v)*p*dx 

# stabilization of p 
a11 = tau*dot(grad(p), grad(q))*dx


# right hand side 
L0 = dot(v, f)*dx 
#FIXME: does not remember stabilization terms on rhs right now
L1 = tau*q*dx  

#backend = EpetraFactory.instance()
backend = PETScFactory.instance()

# Solution vectors 
x0 = backend.createVector()
x1 = backend.createVector()

# Functions
U = Function(V, mesh, x0)
P = Function(Q, mesh, x1)

# create matrices
A00 = assemble(a00, mesh, backend=backend)
A10 = assemble(a10, mesh, backend=backend)
A01 = assemble(a01, mesh, backend=backend)
A11 = assemble(a11, mesh, backend=backend)

# create right hand sides
b0 = assemble(L0, mesh, backend=backend)
b1 = assemble(L1, mesh, backend=backend)

# create block matrix
AA = BlockMatrix(2,2)
AA[0,0] = A00 
AA[1,0] = A10 
AA[0,1] = A01 
AA[1,1] = A11 

file = File("AP00.m"); file << A00
file = File("AP10.m"); file << A10
file = File("AP01.m"); file << A01
file = File("AP11.m"); file << A11

# create block vector 
bb = BlockVector(2)  
bb[0] = b0
bb[1] = b1


# create solution vector
xx = BlockVector(2)  
xx[0] = x0
xx[1] = x1


# create matrices needed in the preconditioner 
c11 = p*q*dx
C11 = assemble(c11, mesh, backend=backend)

# create block preconditioner
BB = SaddlePrec(AA[0,0], C11)


#regularization_parameter = 1/500.0
#xx = precRichardson(BB, AA, xx, bb, regularization_parameter, 10e-8, True, 20)
xx = precondMinRes(BB, AA, xx, bb, 10e-8, True, 200)


