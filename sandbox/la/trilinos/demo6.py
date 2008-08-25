from PyTrilinos import Epetra, AztecOO, TriUtils, ML 
from sys import path 
path.append("../poisson")


import math 
from dolfin import *
#from BlockLinearAlgebra import *
from Krylov import *
from MinRes import *
dolfin_set("linear algebra backend", "Epetra")

class MLPreconditioner: 
    def __init__(self, A): 
        # Sets up the parameters for ML using a python dictionary
        MLList = {
              "max levels"        : 3, 
              "output"            : 0, # as little output as possible
              "smoother: type"    : "ML symmetric Gauss-Seidel",
              "aggregation: type" : "Uncoupled",
              "ML validate parameter list" : False
        }
        ml_prec = ML.MultiLevelPreconditioner(dolfin.down_cast_epetra_matrix(A.instance()).mat(), False)
        ml_prec.SetParameterList(MLList)
        ml_prec.ComputePreconditioner()
        self.ml_prec = ml_prec
    def __mul__(self, b):
        x = b.copy()
        x.zero()
        err = self.ml_prec.ApplyInverse(dolfin.down_cast_epetra_vector(b.instance()).vec(),
                                        dolfin.down_cast_epetra_vector(x.instance()).vec())
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
        if (not b.size() == self.n):
            raise ValueError, "dim(BlockVector) != dim(Prec) "  

        b0 = b[0]
        b1 = b[1]

        x0 = self.M_prec*b0
        x1 = self.N_prec*b1

        xx = BlockVector(2)
        xx[0] = x0 
        xx[1] = x1
        return xx 

class IdentityPrec: 
    def __init__(self):
        pass
    def __mul__(self, b):
        x = b.copy()
        return x 


# Function for no-slip boundary condition for velocity
class Noslip(Function):
    def __init__(self, mesh):
        Function.__init__(self, mesh)

    def eval(self, values, x):
        values[0] = 0.0
        values[1] = 0.0

    def rank(self):
        return 1

    def dim(self, i):
        return 2

# Function for inflow boundary condition for velocity
class Inflow(Function):
    def __init__(self, mesh):
        Function.__init__(self, mesh)

    def eval(self, values, x):
        values[0] = -1.0
        values[1] = 0.0

    def rank(self):
        return 1

    def dim(self, i):
        return 2

# Function for inflow boundary condition for velocity
class Source(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)

    def eval(self, values, x):
        values[0] = x[1]
        values[1] = x[0] 

    def rank(self):
        return 1

    def dim(self, i):
        return 2



mesh = UnitSquare(40,40)
#mesh = Mesh("../../../data/meshes/dolfin-2.xml.gz")
#sub_domains = MeshFunction("uint", mesh, "../../../demo/pde/stokes/taylor-hood/subdomains.xml.gz")

# create Taylor-Hood elements
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
a00 = dot(v,u)*dx + dot(grad(v), grad(u))*dx

# Divergence constraint 
a10 = -div(u)*q*dx

# gradient of p 
a01 = -div(v)*p*dx 

# stabilization of p 
a11 = tau*dot(grad(p), grad(q))*dx


# right hand side 
L0 = dot(v, f)*dx 
#FIXME: does not remember stabilization terms on rhs right now
L1 = tau*q*dx  


# Create functions for boundary conditions
noslip = Noslip(mesh)
inflow = Inflow(mesh)

# No-slip boundary condition for velocity
#bc0 = DirichletBC(noslip, sub_domains, 0)

# Inflow boundary condition for velocity
#bc1 = DirichletBC(inflow, sub_domains, 1)

# Collect boundary conditions
#bcs = [bc0, bc1]


# create matrices
A00 = assemble(a00, mesh)
A10 = assemble(a10, mesh)
A01 = assemble(a01, mesh)
A11 = assemble(a11, mesh)

# create right hand sides
b0 = assemble(L0, mesh)
b1 = assemble(L1, mesh)

# apply bc 
#for bc in bcs: 
#    bc.apply(A00, b0, a00)
#    bc.zero(A01, a00)


# ensure that bc are in start vector
x0 = b0.copy()
x1 = b1.copy()

# create block vector 
bb = BlockVector(2, False)  
bb[0] = b0
bb[1] = b1


# create block matrix
AA = BlockMatrix(2,2, False)
AA[0,0] = A00 
AA[1,0] = A10 
AA[0,1] = A01 
AA[1,1] = A11 

#file = File("A00.m"); file << A00
#file = File("A10.m"); file << A10
#file = File("A01.m"); file << A01
#file = File("A11.m"); file << A11


# create solution vector
xx = BlockVector(2, False)  
xx[0] = x0
xx[1] = x1


# create matrices needed in the preconditioner 
c11 = p*q*dx
C11 = assemble(c11, mesh)
#file = File("C11.m"); file << C11


# create block preconditioner
BB = SaddlePrec(AA[0,0], C11)

#BB = IdentityPrec()

xx, i, rho  = precondMinRes(BB, AA, xx, bb, 10e-8, True, 200)
#xx, i, rho  = precondMinRes(BB, AA, xx, bb, 10e-8, True, 200)
#xx, i, rho  = precondBiCGStab(BB, AA, xx, bb, 10e-8, True, 200)

# Functions
U = Function(V, mesh, xx[0])
P = Function(Q, mesh, xx[1])



# Plot solution
#plot(mesh)
plot(U)
plot(P)
interactive()




