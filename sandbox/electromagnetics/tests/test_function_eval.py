# A script to test the evaluation of a vector function from pydoflin

__author__ = "Evan Lezar (evanlezar@gmail.com)"
__date__ = "2008-08-21"
__copyright__ = "Copyright (C) 2008 Evan Lezar"

from dolfin import *
import numpy as N

def test_zero_function(order = 0):
    """
    This function tests to see that a discrete vector function with zero coefficientts evaluates to zero
    at a random point in a unit square mesh 
    """
    mesh = UnitSquare(2,2)
    element = FiniteElement("Nedelec", "triangle", order)
    
    v = TestFunction(element)
    u = TrialFunction(element)
            
    L = dot(v, u)
    # assemble T to get the size of the system
    T = assemble(L*dx, mesh)
    
    num_dofs = T.size(1)
    
    
    # set the coeficients of all the dofs to 0
    # note that the list of coefficients MUST be a dolfin vector
    h_e = Vector(num_dofs)
    
    # initialize the function
    field = Function(element, mesh, h_e)
            
    point = N.random.rand(1,2)
    
    # randomly initialize H so as to ensure that the function value is set
    H = N.random.rand(1,2)
    
    
    # evaluate H at point        
    field.eval(H, point)
            
    # check to see that all entries of H are zero
    assert ( ~H.all() )
    
    print "Passed"

if __name__ == "__main__":
    order = 0
    test_zero_function(order)