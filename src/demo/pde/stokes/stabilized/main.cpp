// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2006-02-09
// Last changed: 2006-02-10

#include <dolfin.h>
#include "Stokes.h"

using namespace dolfin;

int main()
{
  // Boundary condition
  class : public BoundaryCondition
  {
    void eval(BoundaryValue& value, const Point& p, unsigned int i)
    {
      // Pressure boundary condition, zero pressure at one point
      if ( i == 2 && p.x < DOLFIN_EPS && p.y < DOLFIN_EPS )
      {
	value = 0.0;
	return;
      }
      
      // Velocity boundary condition at inflow
      if ( p.x > (1.0 - DOLFIN_EPS) )
      {
	if ( i == 0 )
	  value = -1.0;
	else
	  value = 0.0;
	return;
      }
      
      // Velocity boundary condition at remaining boundary (excluding outflow)
      if ( p.x > DOLFIN_EPS )
	value = 0.0;
    }
  } bc;

  // Set up problem
  Mesh mesh("dolfin-2.xml.gz");
  Function f = 0.0;
  MeshSize h;
  Stokes::BilinearForm a(h);
  Stokes::LinearForm L(f, h);
  PDE pde(a, L, mesh, bc);

  // Compute solution
  Function u;
  Function p;
  pde.solve(u, p);

  // Save solution to file
  File ufile("velocity.pvd");
  File pfile("pressure.pvd");
  ufile << u;
  pfile << p;
}
