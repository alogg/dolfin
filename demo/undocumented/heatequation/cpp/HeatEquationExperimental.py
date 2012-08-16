# Copyright (C) 2012 Marie E. Rognes (meg@simula.no)
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2012-08-15
# Last changed: 2012-08-16
#
# The bilinear form a(u, v) and linear form L(v) for a cG(1)
# discretization of the time-dependent heat equation.
#
# This is an experimental approximation of HeatEquation.ufl aiming at
# convergence to the latter at which point it should be removed.
#
# Run with python HeatEquationExperimental.py

from ufl import *
from ffc import *

# Define dt
# FIXME: perhaps add to global UFL namespace
dt = Measure("cell")

# Define Dt
# FIXME: perhaps add to global UFL namespace
def Dt(u):
    return u.dx(0)

# Define tensor_product
# FIXME: Add to UFL
from ufl import FiniteElementBase
def tensor_product(*args):
    if not all(isinstance(arg, FiniteElementBase) for arg in args):
        raise RuntimeError, "Expecting arguments of tensor_product operator to be finite elements."
    return TensorProductElement(*args)

# Trial and test spaces for space discretization
Uh = FiniteElement("Lagrange", triangle, 1)
Vh = FiniteElement("Lagrange", triangle, 1)

# Trial and test spaces for time discretization
# FIXME: Awaiting temporal 1D basis function generator
Uk = FiniteElement("CG", interval, 1)
Vk = FiniteElement("DG", interval, 0)

# Trial and test spaces for space-time discretization
# FIXME: Consider short-hand for this. Marie is not convinced about **
U = tensor_product(Uk, Uh)
V = tensor_product(Vk, Vh)

# Trial and test functions
u = TrialFunction(U)
v = TestFunction(V)

# Heat conductivity and source term
kappa = Coefficient(U)
f = Coefficient(U)

# Define product measure
# FIXME: Multiplication of form with measure should be allowed
dxdt = dx*dt

# Bilinear and linear forms
# FIXME: More handlers need to be implemented for Kronecker transform
#a = Dt(u)*v*dxdt + kappa*dot(grad(u), grad(v))*dxdt
a = Dt(u)*v*dxdt + sum(kappa*u.dx(i)*v.dx(i)
                       for i in range(1, Uh.cell().d))*dxdt
L = f*v*dxdt

parameters = default_parameters()
parameters["format"] = "dolfin"

# Compile forms
compile_form_kronecker(a, prefix="Heat_2D", parameters=parameters)
compile_form_kronecker(L, prefix="Heat_2D_rhs", parameters=parameters)
