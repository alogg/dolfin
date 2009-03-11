from dolfin import *

def test1():
    "Test time-stepping with initial data expression."

    # Initializations
    mesh = UnitSquare(32, 32)
    V = FunctionSpace(mesh, "Lagrange", 1)
    u0 = Function(V, "sin(x[0])")
    u1 = Function(V)

    # Time stepping
    for i in range(10):

        print i

        # Solve for u1
        u1.vector()

        # Assign u0 = u1 (works fine)
        u0.assign(u1)

def test2():
    "Test time-stepping with time-dependent coefficient expression."

    # Initializations
    mesh = UnitSquare(32, 32)
    V = FunctionSpace(mesh, "Lagrange", 1)
    w0 = Function(V)
    w1 = Function(V, "sin(t*x[0])")

    # Time stepping
    for i in range(10):

        print i

        # Update w1
        w1.t = float(i)
        
        # Solve for u

        # One of these needed

        # Assign w0 = w1 (does not work)
        w0.assign(w1)
        #w0 = interpolate(w1, V)
        #w0 = project(w1, V)

test1()
test2()
