// Copyright (C) 2005 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#include <dolfin/dolfin_log.h>
#include <dolfin/dolfin_math.h>
#include <dolfin/Alloc.h>
#include <dolfin/ODE.h>
#include <dolfin/NewTimeSlab.h>
#include <dolfin/NewMethod.h>
#include <dolfin/NewtonSolver.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
NewtonSolver::NewtonSolver(ODE& ode, NewTimeSlab& timeslab, const NewMethod& method)
  : TimeSlabSolver(ode, timeslab, method), f(0), A(ode, timeslab, method)
{
  // Initialize local array
  f = new real[method.qsize()];
}
//-----------------------------------------------------------------------------
NewtonSolver::~NewtonSolver()
{
  // Delete local array
  if ( f ) delete [] f;
}
//-----------------------------------------------------------------------------
void NewtonSolver::start()
{
  // Get size of system
  int nj = static_cast<int>(ts.nj);

  // Initialize increment vector
  dx.init(nj);

  // Initialize right-hand side
  b.init(nj);

  // Initialize Jacobian matrix
  A.init(dx, dx);

  // Compute Jacobian
  A.update(ts.starttime());

  A.disp();
}
//-----------------------------------------------------------------------------
real NewtonSolver::iteration()
{
  // Evaluate b = -F(x) at current x
  beval();
  
  // Solve linear system F for dx
  solver.solve(A, dx, b);
   
  // Get array containing the increments (assumes uniprocessor case)
  real* dxvals = dx.array();

  /*
  real* bvals = b.array();
  for (uint j = 0; j < ts.nj; j++)
    dxvals[j] = bvals[j];
  b.restore(bvals);
  */

  // Update solution x -> x - dx
  for (uint j = 0; j < ts.nj; j++)
    ts.jx[j] += dxvals[j];

  // Compute maximum increment
  real max_increment = 0.0;
  for (uint j = 0; j < ts.nj; j++)
  {
    const real increment = fabs(dxvals[j]);
    if ( increment > max_increment )
      max_increment = increment;
  }

  // Restore array
  dx.restore(dxvals);

  return max_increment;
}
//-----------------------------------------------------------------------------
void NewtonSolver::beval()
{
  // Get array of values for F (assumes uniprocessor case)
  real* bvals = b.array();

  // Reset dof
  uint j = 0;

  // Reset current sub slab
  int s = -1;

  // Reset elast
  ts.elast = -1;

  // Iterate over all elements
  for (uint e = 0; e < ts.ne; e++)
  {
    // Cover all elements in current sub slab
    s = ts.cover(s, e);

    // Get element data
    const uint i = ts.ei[e];
    const real a = ts.sa[s];
    const real b = ts.sb[s];
    const real k = b - a;

    // Get initial value for element
    const int ep = ts.ee[e];
    const real x0 = ( ep != -1 ? ts.jx[ep*method.nsize() + method.nsize() - 1] : ts.u0[i] );

    // Evaluate right-hand side at quadrature points of element
    ts.feval(f, s, e, i, a, b, k);
    //cout << "f = "; Alloc::disp(f, method.qsize());

    // Update values on element using fixed point iteration
    method.update(x0, f, k, bvals + j);
    
    // Subtract current values
    for (uint n = 0; n < method.nsize(); n++)
      bvals[j + n] -= ts.jx[j + n];

    // Update dof
    j += method.nsize();
  }

  // Restor array
  b.restore(bvals);
}
//-----------------------------------------------------------------------------
