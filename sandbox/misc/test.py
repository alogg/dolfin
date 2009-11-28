from dolfin import *
from numpy  import array, zeros

f1 = Expression("2", degree=2)
print "-------------------"
f2,f3 = Expressions("2", "3", element=FiniteElement("CG", tetrahedron, 2))
#f2,f3 = Expressions("2", "3", degree=2)


#u20 = zeros(1)
#print "Eval"

print "Test A ",  f1(0.31, 0.32, 0.33)
print "Test B ",  f2(0.31, 0.32, 0.33)
#print "End Eval"
#f2(0.31, 0.32, 0.33, values = u20)

