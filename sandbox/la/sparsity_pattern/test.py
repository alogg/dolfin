"""A simple test for displaying the sparsity pattern for uBLAS and MTL4 backends"""

__author__ = "Johan Hake (hake@simula.no)"
__date__ = "2009-02-18"
__copyright__ = "Copyright (C) 2008 Johan Hake"
__license__  = "GNU LGPL Version 2.1"

from dolfin import *

dolfin_set("linear algebra backend","uBLAS")

V = FunctionSpace(UnitSquare(3,3),"Lagrange", 1)
a = dot(grad(TestFunction(V)),grad(TrialFunction(V)))*dx

A = assemble(a)

cols, rows, values = A.data()

print "Sparsity pattern with uBLAS backend:"
print "cols",cols
print "rows",rows
print "values",values
print "nnz",len(values)

dolfin_set("linear algebra backend","MTL4")

A = assemble(a)

cols, rows, values = A.data()

print "\nSparsity pattern with MTL4 backend:"
print "cols",cols
print "rows",rows
print "values",values
print "nnz",len(values)

