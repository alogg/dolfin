__author__ = "Ilmar Wilbers (ilmarw@simula.no)"
__date__ = "2008-06-04"
__copyright__ = "Copyright (C) 2007 Ilmar Wilbers"
__license__  = "GNU LGPL Version 2.1"

from dolfin import *
from time import time
import sys

def make_form(name):
    execfile(name + ".form", globals())
    try:
        return a
    except:
        print "No object 'a' to return in file %s.form" %name
        return None

def bench_form(form, mesh, nr_reasm=1, coeffs=[]):
    totaltime = 0.0
    for i in range(nr_reasm):
        t0 = time()
        A = assemble(form, mesh, coefficients=coeffs)
        totaltime = time() - t0
    return totaltime/float(nr_reasm)

def make_mesh(name, dim=2):
    if dim == 3:
        N = 32
        mesh = UnitCube(N, N, N)
        return mesh
    else: # Assume 3D
        N = 256
        mesh = UnitSquare(N, N)
        return mesh

def print_table(values, title):
    "Print nicely formatted table"
    m = max([key[0] for key in values]) + 2
    n = max([key[1] for key in values]) + 2

    table = []
    for i in range(m):
        table.append(["" for j in range(n)])

    for i in range(m - 1):
        table[i + 1][0] = str(values[(i, 0)][0]).split(" ")[0]

    for j in range(n - 1):
        table[0][j + 1] = str(values[(0, j)][1]).split(" ")[0]

    for i in range(m - 1):
        for j in range(n - 1):
            table[i + 1][j + 1] = "%.5g" % values[(i, j)][2]

    table[0][0] = title

    column_sizes = [max([len(table[i][j]) for i in range(m)]) for j in
    range(n)]
    row_size = sum(column_sizes) + 3*(len(column_sizes) - 1) + 2

    print ""
    for i in range(m):
        print " " + "-"*row_size
        print "|",
        for j in range(n):
            print table[i][j] + " "*(column_sizes[j] - len(table[i][j])),
            print "|",
        print ""
    print " " + "-"*row_size
    print ""

if __name__ == "__main__":
    try:
        num_reasm = int(sys.argv[1])
    except:
        print 'Usage: %s number_of_repetitions [1]' %sys.argv[0]
        num_reasm = 1

dolfin_set("output destination", "silent")
backends = ["uBLAS", "PETSc", "Epetra", "Assembly"]
# forms consist of name of the form (same as file)
# and a list of arguments to assembler (given to the constructor in cpp)
# as a tuple.
forms = [("Elasticity3D", None),
         ("PoissonP1", None),
         ("PoissonP2", None),
         ("PoissonP3", None),
         ("THStokes2D", None),
#         ("NSEMomentum3D", [1., 1., 1., 1., 1.]),
#         ("StabStokes2D", [1.]),
         ]
results = {}
for i in range(len(backends)):
    backend = backends[i]
    dolfin_set("linear algebra backend", backends[i])
    for j in range(len(forms)):
        form = forms[j]
        dim = 2 if not form[0].find("3D") > -1 else 3
        a = make_form(form[0])
        m = make_mesh(form[0], dim)
        print "Assembling %s with %s" %(form[0], backend)
        results[(i, j)] = (backend, form[0], bench_form(a, m,
                                                        nr_reasm=num_reasm,
                                                        coeffs=form[1]))

print_table(results, "Backend/Form")
