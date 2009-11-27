"""This demo shows the intersection of the boundary of a unit square
(omega1) with a unit circle (omega2) rotating around the center of the
square."""

__author__ = "Andre Massing (masssing@simula.no)"
__date__ = ""
__copyright__ = "Copyright (C) 2009 Andre Massing"
__license__  = "GNU LGPL Version 2.1"

from dolfin import *
from numpy import *

#modify to check for CGAL
#if not has_gts():
#    print "DOLFIN must be compiled with GTS to run this demo."
#    exit(0)

# Create meshes (omega0 overlapped by omega1)
#mesh0 = UnitSquare(3, 3)
#mesh1 = UnitCircle(5)
mesh0 = UnitCube(3, 3, 2)
mesh1 = UnitSphere(3)

print "Cube"
print "Cells:"
#print mesh0.cells() 
#print mesh0.coordinates()
print mesh0.coordinates()[mesh0.cells()[102]]

print 

#Create intersection operator for mesh0
iop = IntersectionOperator(mesh0)

# Access mesh geometry
x = mesh1.coordinates()

# Move and scale circle
x *= 0.5
x += 1.0

print "Sphere"
print "Cells:"
print mesh1.coordinates()[mesh1.cells()[79]]
print mesh1.cells() 
print mesh1.coordinates()

cells = iop.all_intersected_entities(mesh1)
