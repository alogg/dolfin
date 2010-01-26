// Copyright (C) 2009 Andre Massing
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2010-10-08
// Last changed: 2010-01-14

#include <dolfin.h>
#include <math.h>

using namespace dolfin;

#ifdef HAS_CGAL

int main()
{
  //Mesh omega2 overlaps omega1.
  UnitCube omega1(10, 10, 10);
  UnitSphere omega2(10);

  // Access mesh geometry
  MeshGeometry& geometry = omega2.geometry();

  // Move and scale circle
  for (VertexIterator vertex(sphere); !vertex.end(); ++vertex)
  {
    double* x = geometry.x(vertex->index());
    x[0] = 0.5*x[0] + 1.0;
    x[1] = 0.5*x[1] + 1.0;
  }

  //Compute cells of omega1 which are partially intersected by omega2 by
  //intersecting boundary of omega2 with omega1.
  BoundaryMesh boundary(omega2);

  // Compute intersection with boundary of square
  //typedef for std::set<unsigned int>
  uint_set cells;
  omega0.all_intersected_entities(boundary, cells);

}

#else

int main()
{
  info("DOLFIN must be compiled with CGAL to run this demo.");
  return 0;
}

#endif
