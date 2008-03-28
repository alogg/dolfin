// Copyright (C) 2008 Solveig Bruvoll and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-03-28
// Last changed: 2007-03-28
//
// Transfinite mean value interpolation (ALE mesh smoothing)

#include <dolfin.h>
#include <string.h>
#include <cmath>

using namespace dolfin;

real dist(const real* x, const real* y, unsigned int dim)
{
  real s = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    s += (x[i] - y[i])*(x[i] - y[i]);
  return sqrt(s);
}

void meanValueH(const real* old_x, real* new_x, unsigned int dim, Mesh& new_boundary,
                Mesh& mesh, MeshFunction<unsigned int>& vertex_map)
{
  cout << old_x[0] << " " << old_x[1] << endl;

  const unsigned int num_vertices = new_boundary.topology().dim() + 1;
  real** old_p = new real * [num_vertices];
  real** new_p = new real * [num_vertices];

  for (CellIterator c(new_boundary); !c.end(); ++c)
  {
    for (VertexIterator v(*c); !v.end(); ++v)
    {
      old_p[v.pos()] = mesh.geometry().x(vertex_map(*v));
      new_p[v.pos()] = v->x();
    }

    //old_p[0][0]
    //new_p

  }

  // for (unsigned int i = 0; i < dim; i++)
  // new_x[i] = ...

  delete [] old_p;
  delete [] new_p;
}

void deform(Mesh& mesh, Mesh& new_boundary, MeshFunction<unsigned int>& vertex_map, MeshFunction<unsigned int>& cell_map)
{
  // Extract old coordinates
  const unsigned int dim = mesh.geometry().dim();
  const unsigned int size = mesh.numVertices()*dim;
  real* old_x = new real[size];
  real* new_x = mesh.geometry().x();
  memcpy(old_x, new_x, size*sizeof(real));

  // Iterate over coordinates in mesh
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    meanValueH(old_x + v->index()*dim, new_x + v->index()*dim, dim, new_boundary,
               mesh, vertex_map);
  }

  delete [] old_x;
}

int main()
{
  UnitCube mesh(3, 3, 3);
  plot(mesh);

  MeshFunction<unsigned int> vertex_map;
  MeshFunction<unsigned int> cell_map;
  BoundaryMesh boundary(mesh, vertex_map, cell_map);

  // Deform boundary
  for (VertexIterator v(boundary); !v.end(); ++v)
  {
    real* x = v->x();
    x[0] = x[0] + x[2];
  }
  plot(boundary);

  // Deform mesh
  deform(mesh, boundary, vertex_map, cell_map);
  plot(mesh);

  return 0;
}
