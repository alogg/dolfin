// Copyright (C) 2006 Garth N. Wells.
// Licensed under the GNU GPL Version 2.
//
// First added:  2006-06-01
// Last changed: 2006-06-06

#include <dolfin/dolfin_log.h>
#include <dolfin/uBlasLUSolver.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
uBlasLUSolver::uBlasLUSolver() : LinearSolver()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
uBlasLUSolver::~uBlasLUSolver()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
dolfin::uint uBlasLUSolver::solve(const uBlasSparseMatrix& A, DenseVector& x, 
    const DenseVector& b)
{
  dolfin_warning("LU solver will be used. This may be slow and consume a lot of memory.");
        
  // FIXME: implement renumbering scheme to speed up LU solve
  A.solve(x, b);
  return 1;
}
//-----------------------------------------------------------------------------
