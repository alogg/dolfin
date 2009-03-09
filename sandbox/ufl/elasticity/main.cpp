// Copyright (C) 2006-2009 Johan Jansson and Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2008
//
// First added:  2006-02-07
// Last changed: 2009-02-25
//
// This demo program solves the equations of static
// linear elasticity for a gear clamped at two of its
// ends and twisted 30 degrees.

#include <dolfin.h>
#include "Elasticity.h"
#include "ElasticityUFL.h"

using namespace dolfin;

int main()
{
  // Read mesh and create function space
  //Mesh mesh("../../../data/meshes/gear.xml.gz");
  UnitCube mesh(30, 30, 30);
  ElasticityFunctionSpace V(mesh);
  ElasticityUFLFunctionSpace V_ufl(mesh);

  // Set elasticity parameters
  double E  = 10.0;
  double nu = 0.3;
  Constant mu(E / (2*(1 + nu)));
  Constant lambda(E*nu / ((1 + nu)*(1 - 2*nu)));

  // Create right-hand side
  Constant f(3, 0.0);

  // Set up forms
  ElasticityBilinearForm a(V, V);
  a.mu = mu; a.lmbda = lambda;
  ElasticityLinearForm L(V);
  L.f = f;

  ElasticityUFLBilinearForm a_ufl(V_ufl, V_ufl);
  a_ufl.mu = mu; a_ufl.lmbda = lambda;
  ElasticityLinearForm L_ufl(V_ufl);
  L_ufl.f = f;

  Matrix A, A_ufl;
  tic();
  assemble(A, a);
  double time0 = toc();

  tic();
  assemble(A_ufl, a_ufl);
  double time1 = toc();

  cout << "Assembly timings (old FFC, new UFL) " << time0 << "  " << time1 << endl;

  return 0;
}
