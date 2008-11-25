from ffc import *

# Reserved variables for forms
(a, L, M) = (None, None, None)

# Reserved variable for element
element = None

# Copyright (c) 2005 Johan Hoffman 
# Licensed under the GNU GPL Version 2
#
# Modified by Anders Logg 2006
#
# First added:  2005
# Last changed: 2006-03-28
#
# The momentum equation for the incompressible 
# Navier-Stokes equations using cG(1)cG(1)
#
# Compile this form with FFC: ffc NSEMomentum2D.form.

name = "NSEMomentum3D"
scalar = FiniteElement("Lagrange", "tetrahedron", 1)
vector = VectorElement("Lagrange", "tetrahedron", 1)
constant_scalar = FiniteElement("Discontinuous Lagrange", "tetrahedron", 0)
constant_vector = VectorElement("Discontinuous Lagrange", "tetrahedron", 0)

v      = TestFunction(vector)      # test basis function
U      = TrialFunction(vector)     # trial basis function
#uc     = Function(vector) # linearized velocity
um     = Function(constant_vector) # cell mean linearized velocity
u0     = Function(vector)          # velocity from previous time step
f      = Function(vector)          # force term
p      = Function(scalar)          # pressure
delta1 = Function(constant_scalar) # stabilization parameter
delta2 = Function(constant_scalar) # stabilization parameter

k  = Constant("tetrahedron") # time step
nu = Constant("tetrahedron") # viscosity

#uc = mean(uc)   # cell mean value of linearized velocity
#uc = uc   # cell mean value of linearized velocity

i0 = Index()    # index for tensor notation
i1 = Index()    # index for tensor notation
i2 = Index()    # index for tensor notation

# Galerkin discretization of bilinear form  
G_a  = (dot(v, U) + k*nu*0.5*dot(grad(v), grad(U)) + 0.5*k*v[i0]*um[i1]*U[i0].dx(i1))*dx
# Least squares stabilization of bilinear form  
SD_a = (delta1*k*0.5*um[i1]*v[i0].dx(i1)*um[i2]*U[i0].dx(i2) + delta2*k*0.5*div(v)*div(U))*dx

# Galerkin discretization of linear form  
G_L  = (dot(v, u0) + k*dot(v, f) + k*div(v)*p - k*nu*0.5*dot(grad(v), grad(u0)) - 0.5*k*v[i0]*um[i1]*u0[i0].dx(i1))*dx
# Least squares stabilization of linear form
SD_L = (- delta1*k*0.5*um[i1]*v[i0].dx(i1)*um[i2]*u0[i0].dx(i2) - delta2*k*0.5*div(v)*div(u0))*dx

# Bilinear and linear forms
a = G_a + SD_a
L = G_L + SD_L


compile([a, L, M, element], "NSEMomentum3D", options={'language': 'dolfin', 'blas': False, 'form_postfix': True, 'precision': '15', 'cache_dir': None, 'cpp optimize': False, 'split_implementation': False, 'quadrature_points': False, 'output_dir': '.', 'representation': 'tensor', 'shared_ptr': False, 'optimize': False}, global_variables=globals())
