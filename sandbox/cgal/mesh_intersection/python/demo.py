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

print "Sphere"
print "Cells:"
print mesh1.coordinates()[mesh1.cells()[79]]
#print mesh1.cells() 
#print mesh1.coordinates()

#Create intersection operator for mesh0
iop = IntersectionOperator(mesh0)

# Access mesh geometry
x = mesh1.coordinates()

# Move and scale circle
x *= 0.5
x += 1.0

# Iterate over angle
theta = 0.0
dtheta = 0.01*DOLFIN_PI
intersection = MeshFunction("uint", mesh0, mesh0.topology().dim())
_first = False

cells = iop.all_intersected_entities(mesh1)
p = plot(intersection, rescale=False)
p.add_polygon([[0, 0, -0.01], [1, 0, -0.01], [1, 1, -0.01], [0, 1, -0.01], [0, 0, -0.01]])
p.ren.ResetCamera()

while theta >  2*DOLFIN_PI:

    # Compute intersection with boundary of square
    cells = iop.all_intersected_entities(mesh1)
#    boundary =BoundaryMesh(mesh1)
#    cells = iop.all_intersected_entities(boundary)

    # Mark intersected values
    intersection.values()[:] = 0
    intersection.values()[cells] = 1

    # Plot intersection
    if _first:
        p = plot(intersection, rescale=False)
        p.add_polygon([[0, 0, -0.01], [1, 0, -0.01], [1, 1, -0.01], [0, 1, -0.01], [0, 0, -0.01]])
        p.ren.ResetCamera()
        _first = False
    else:
        plot(intersection)

#    interactive()


    # Rotate circle around (0.5, 0.5)
    xr = x[:, 0] - 0.5
    yr = x[:, 1] - 0.5
    x[:,0] = 0.5 + (cos(dtheta)*xr - sin(dtheta)*yr)
    x[:,1] = 0.5 + (sin(dtheta)*xr + cos(dtheta)*yr)

    theta += dtheta

# Hold plot
interactive()
