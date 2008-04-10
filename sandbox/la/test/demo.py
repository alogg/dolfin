
from dolfin import *
from sys import path 
path.append("../poisson")


from BlockLinearAlgebra import *

# Create mesh and finite element
N = 4 
mesh = UnitSquare(N,N)
mesh.disp()
CG = FiniteElement("Lagrange", "triangle", 1)
DG = FiniteElement("DG", "triangle", 0)

# Source term
class Source1(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)
    def eval(self, values, x):
        values[0] = exp(x[0]*x[1]) 

# Source term
class Source2(Function):
    def __init__(self, element, mesh):
        Function.__init__(self, element, mesh)
    def eval(self, values, x):
        values[0] = cos(x[0]*x[1]) 


v = TestFunction(CG)
u = TrialFunction(DG)
w = TestFunction(DG)
f1 = Source1(CG, mesh)
f2 = Source2(DG, mesh)
f3 = Source2(CG, mesh)

a = u*v*dx + dot(grad(v), grad(u))*dx 
L1 = v*f1*dx 
L2 = w*f2*dx 
L3 = v*f2*dx 

# Assemble matrices
A = assemble(a, mesh)
b1 = assemble(L1, mesh)
b2 = assemble(L2, mesh)
b3 = assemble(L3, mesh)

file = File("A.m"); file << A; 
file = File("b1.m"); file << b1;  
file = File("b2.m"); file << b2;  
file = File("b3.m"); file << b3;  

print ""
print "testing Matrix Vector Product"
print " testing Vector (transposed)", 
x = b2.copy()
x.zero()
A.mult(b1, x, True)
print x.inner(x)

x.zero()
print " testing BlockVector (transposed)", 
AA = BlockMatrix(1,1); AA[0,0] = A; 
bb1 = BlockVector(1) ; bb1[0] = b1;  
xx = BlockVector(1); xx[0] = x; 

AA.prod(bb1, xx, True);
print xx.inner(xx); 


print ""
print "testing Matrix Vector Product"
print " testing Matrix::mult ", 
x = b1.copy()
x.zero()
A.mult(b2, x, False)
print x.inner(x)
print " testing Matrix::__mul__ ", 
x = A*b2
print x.inner(x)


print " testing BlockVector", 

x.zero()
bb2 = BlockVector(1) ; bb2[0] = b2;  
xx = BlockVector(1); xx[0] = x; 

AA.prod(bb2, xx, False);
print xx.inner(xx); 


print ""
print "testing Vector Addition "
tmp = b1.copy()
print " testing Vector", 
tmp += b3 
print tmp.inner(tmp)

tmp = b1.copy()
bb1 = BlockVector(1) ; bb1[0] = tmp;  
bb3 = BlockVector(1) ; bb3[0] = b3;  
bb1 += bb3
print " testing BlockVector", bb1.inner(bb1)

print ""
print "testing Vector Subtraction "
tmp = b1.copy()
print " testing Vector", 
tmp -= b3 
print tmp.inner(tmp)

tmp = b1.copy()
bb1 = BlockVector(1) ; bb1[0] = tmp;  
bb3 = BlockVector(1) ; bb3[0] = b3;  
bb1 -= bb3
print " testing BlockVector", bb1.inner(bb1)


print ""
print "testing Vector Multiplication "
tmp = b1.copy()
print " testing Vector", 
tmp *=3.7 
print tmp.inner(tmp)

tmp = b1.copy()
bb1 = BlockVector(1) ; bb1[0] = tmp;  
bb3 = BlockVector(1) ; bb3[0] = b3;  
bb1 *= 3.7
print " testing BlockVector", bb1.inner(bb1)












