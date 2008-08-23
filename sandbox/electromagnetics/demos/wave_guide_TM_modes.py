""" 
This demo demonstrates the calculation and visualization of a TM (Transverse Magnetic) cutoff mode of a rectangular waveguide.
Note that for simplicity, the relative permiability (\mu_r) and permitivity (\epsilon_r) are taken as unity.

For more information regarding waveguides see 

http://www.ee.bilkent.edu.tr/~microwave/programs/magnetic/rect/info.htm

"""
__author__ = "Evan Lezar evanlezar@gmail.com"
__date__ = "2008-08-22"
__copyright__ = "Copyright (C) 2008 Evan Lezar"


from dolfin import *
import numpy as N

# specify the waveguide width (a) and height (b) in metres
mesh_a = 1
mesh_b = 0.5

# specify the mode of interest. moi = 1 : dominant (TM_{11}) mode
moi = 1

# create the mesh using a Rectangle
nx = int(mesh_a/mesh_b)
if nx == 0:
    nx = 1
            
mesh = Rectangle(0, mesh_a, 0, mesh_b, nx, 1, 0)

# refine if desired
mesh.refine()
mesh.refine()

# define the finite element.  For vector electromagnetic problems Nedelec vector elements are used.
# Specify the degree of the approximation
degree = 0
element = FiniteElement("Nedelec", "triangle", degree)

# define the test and trial functions
v = TestFunction(element)
u = TrialFunction(element)

# define the forms - gererates an generalized eigenproblem of the form 
# [S]{h} = k_o^2[T]{h}
# with the eigenvalues k_o^2 representing the square of the cutoff wavenumber and the corresponding right-eigenvector giving 
# the coefficients of the discrete system used to obtain the approximate field anywhere in the domain   
 
a = dot(curl_t(v), curl_t(u))*dx
L = dot(v, u)*dx

# Assemble the system matrices
S = assemble(a, mesh)
T = assemble(L, mesh)

# now solve the eigen system
# TODO: Use a sparse iterative solver - scipy 0.7 includes ARPACK

# as a first approximation rewrite the system as [A]{h} = k_o^2{h} with [A] = inv([T])*[S]
A = N.dot(N.linalg.inv(T.array()), S.array())

# U is a vector of eigenvalues and V is a matrix of corresponding eigenvectors
(U, V) = N.linalg.eig(A)

# the result should have real eigenvalues but due to rounding errors, some of the resultant eigenvalues are very small complex values.
# only consider the real part
U = U.real

# Now, the system contains a number of zero eigenvalues (near zero due to rounding) which are eigenvalues corresponding to the null-space
# of the curl operator and are a mathematical construct and do not represent physically realizable modes.  These are called spurious modes.
# So, we need to identify the smallest, non-zero eigenvalue of the system - which corresponds with cutoff wavenumber of the the dominant 
# cutoff mode.

# we need to be able to identify the corresponding eigenvector
sorted_index = U.argsort()
first_mode_index = N.where(U[sorted_index] > 1e-3)[0][0]

# get the square of the cutoff wavenumber and the basis function coefficients
k_o_squared = U[sorted_index[first_mode_index+moi-1]]
h_e_numpy = V[:,sorted_index[first_mode_index+moi-1]].real

print "Cutoff wavenumber squared: %f" % k_o_squared, 

# now to visualize the magnetic field we need to calculate the field at a number of points
# first define a discrete function using the eigenvector values as basis function coefficients
# NOTE:  The coefficients need to be passed to the Function constructor as a dolfin Vector

h_e_dolfin = Vector(len(h_e_numpy))
for k in range(len(h_e_numpy)):
    h_e_dolfin[k] = h_e_numpy[k]
    
# initialize the function
magnetic_field = Function(element, mesh, h_e_dolfin)

# now specify the points where the field must be calculated and calculate
# number of points per unit length
n = 20 

# allocate numpy arrays for the magnetic field (H - x and y components) and the position where the field is calculated 
H = N.zeros(((n+1)*(n+1), 2),dtype=N.float64)
XY = N.zeros(((n+1)*(n+1), 2),dtype=N.float64)
for i in range(n+1):
    for j in range(n+1):
        p_idx = i*(n+1) + j # the index of the point in the array
        XY[p_idx, 0] = float(mesh_a*i)/n
        XY[p_idx, 1] = float(mesh_b*j)/n
        
        # evaluate the magnetic field.  Result is stored in H[p_idx,:]
        magnetic_field.eval(H[p_idx,:], XY[p_idx,:])

# now plot the result
import pylab as P
P.figure()
P.quiver(XY[:,0],XY[:,1],H[:,0],H[:,1])


P.xlabel('x')
P.ylabel('y')

P.axis('equal')

P.show()