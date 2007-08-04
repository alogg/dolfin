from dolfin import *
from ffc import *

# Rename wrapped assemble function
#assemble_cpp = assemble

def foo_assemble(form, mesh):
    print "Assembling..."

    # Compile form
    compiled_form = jit(form, language="dolfin")

    # Assemble compiled form
    rank = compiled_form.rank()
    if rank == 0:
        return assemble_cpp(form, mesh)
    elif rank == 1:
        b = Vector()
        assemble_cpp(b, form, mesh)
        return b
    elif rank == 2:
        A = Matrix()
        assemble_cpp(A, form, mesh)
        return A
    else:
        raise RuntimeError, "Unable to assemble tensors of rank %d." % rank

element = FiniteElement("Lagrange", "triangle", 1)

v = TestFunction(element)
u = TrialFunction(element)

a = dot(grad(v), grad(u))*dx

mesh = UnitSquare(10, 10)
#A = assemble(a, mesh)

A = Matrix()

print isinstance(A, Matrix)
print isinstance(A, GenericMatrix)
print isinstance(A, GenericTensor)

compiled_form = jit(a)
#print isinstance(compiled_form, Form)

#print Form

coefficients = ArrayFunctionPtr()
cell_domains = MeshFunction("uint")
exterior_facet_domains = MeshFunction("uint")
interior_facet_domains = MeshFunction("uint")
assemble(A, a, mesh, coefficients, cell_domains, exterior_facet_domains, interior_facet_domains)

#assemble(1, 1, 1, 1)
