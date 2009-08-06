from dolfin import *

def solve(mesh_file, degree):
    "Solve on given mesh file and degree of function space"

    # Create mesh and define function space
    mesh = Mesh(mesh_file);
    V = FunctionSpace(mesh, "CG", degree)

    # Define variational problem
    v = TestFunction(V)
    u = TrialFunction(V)
    f = Function(V, "sin(x[0])")
    a = (grad(v), grad(u)) + (v, u)
    L = (v, f)

    # Compute solution
    problem = VariationalProblem(a, L)
    problem.parameters["linear_solver"] = "iterative"
    problem.parameters["krylov_solver"]["relative_tolerance"] = 1e-20
    u = problem.solve()

    # Return norm of solution vector
    return u.vector().norm()

def print_reference(results):
    "Print nicely formatted values for gluing into code as a reference"
    print "reference = {",
    for (i, result) in enumerate(results):
        if i > 0:
            print "             ",
        print "(\"%s\", %d): %.16g" % result,
        if i < len(results) - 1:
            print ","
        else:
            print "}"

def check_results(results, reference, tol):
    "Compare results with reference"

    if not MPI.process_number() == 0:
        return

    print "Checking results"
    print "----------------"

    for (mesh_file, degree, norm) in results:
        print "(%s, %d):\t" % (mesh_file, degree),
        if (mesh_file, degree) in reference:
            ref = reference[(mesh_file, degree)]
            diff =  abs(norm - ref)
            if diff < tol:
                print "OK",
            else:
                print "*** WRONG",
            print "(norm = %.16g, reference = %.16g, diff = %.16g)" % (norm, ref, diff)
        else:
            print "missing reference"

# Reference values for norm of solution vector
reference = { ("unitsquare.xml.gz", 1): 7.821707395007537 ,
              ("unitsquare.xml.gz", 2): 15.18829494599347 ,
              ("unitsquare.xml.gz", 3): 22.55234140275229 ,
              ("unitcube.xml.gz", 1): 3.647913575216382 ,
              ("unitcube.xml.gz", 2): 8.523874310611367 ,
              ("unitcube.xml.gz", 3): 14.55432230797502 }

# Mesh files and degrees to check
mesh_files = ["unitsquare.xml.gz", "unitcube.xml.gz"]
degrees = [1, 2, 3]

## Iterate over test cases and collect results
results = []
for mesh_file in mesh_files:
    for degree in degrees:
        norm = solve(mesh_file, degree)
        results.append((mesh_file, degree, norm))

# Uncomment to print results for use as reference
#print_reference(results)

# Check results
check_results(results, reference, 1e-10)
