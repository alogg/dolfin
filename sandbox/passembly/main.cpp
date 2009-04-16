// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-07
// Last changed: 2008-12-26
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = 500*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)
//
// and boundary conditions given by
//
//     u(x, y) = 0 for x = 0 or x = 1

#include <dolfin.h>
#include "Poisson.h"
#include "Poisson3D.h"

using namespace dolfin;

// Source term
class Source : public Function
{
  void eval(double* values, const double* x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    values[0] = 500.0*exp(-(dx*dx + dy*dy) / 0.02);
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const double* x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
  }
};

int test2D()
{
  // Create mesh and function space
  Mesh mesh("unitsquare.xml.gz");
  PoissonFunctionSpace V(mesh);

  PoissonBilinearForm a(V, V);
  UFC ufc(a);
  
  UFC_PoissonBilinearForm_dof_map_0 ufc_dof_map;
  DofMap dofmap(ufc_dof_map, mesh);
  dofmap.build(ufc, mesh);

  char filename[100];
  sprintf(filename, "unitsquare_part_%d.xml", dolfin::MPI::process_number());
  File file(filename, true);
  file << mesh;

  return 0;
}

int test3D()
{
  Mesh mesh("unitcube.xml.gz");
  Poisson3DFunctionSpace V(mesh);
  Poisson3DBilinearForm a(V, V);
  UFC ufc(a);
  UFC_Poisson3DBilinearForm_dof_map_0 ufc_dof_map;
  DofMap dofmap(ufc_dof_map, mesh);
  dofmap.build(ufc, mesh);

  char filename[100];
  sprintf(filename, "unitcube_part_%d.xml", dolfin::MPI::process_number());
  File file(filename, true);
  file << mesh;

  return 0;
}

int main()
{
  return test2D();
}
