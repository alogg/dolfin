// Copyright (C) 2004-2005 Johan Jansson.
// Licensed under the GNU GPL Version 2.
//
// Modified by Anders Logg 2005-2006.
// Modified by Johan Hoffman 2005.
// Modified by Andy R. Terrel 2005.
// Modified by Garth N. Wells 2005-2006.
//
// First added:  2005-12-02
// Last changed: 2006-08-15

#ifndef __PETSC_KRYLOV_SOLVER_H
#define __PETSC_KRYLOV_SOLVER_H

#ifdef HAVE_PETSC_H

#include <dolfin/constants.h>
#include <dolfin/Parametrized.h>
#include <dolfin/Preconditioner.h>
#include <dolfin/KrylovMethod.h>
#include <dolfin/PETScPreconditioner.h>
#include <dolfin/PETScLinearSolver.h>
#include <dolfin/PETScManager.h>

namespace dolfin
{

  /// Forward declarations
  class PETScMatrix;
  class PETScVector;
  class PETScKrylovMatrix;

  /// This class implements Krylov methods for linear systems
  /// of the form Ax = b. It is a wrapper for the Krylov solvers
  /// of PETSc.
  
  class PETScKrylovSolver : public PETScLinearSolver, public Parametrized
  {
  public:

    /// Create Krylov solver with PETSc default method and PETScPreconditioner
    PETScKrylovSolver();

    /// Create Krylov solver for a particular method with default PETSc PETScPreconditioner
    PETScKrylovSolver(KrylovMethod method);

    /// Create Krylov solver with default PETSc method and a particular PETScPreconditioner
    PETScKrylovSolver(Preconditioner pc);

    /// Create Krylov solver with default PETSc method and a particular PETScPreconditioner
    PETScKrylovSolver(PETScPreconditioner& PETScPreconditioner);

    /// Create Krylov solver for a particular method and PETScPreconditioner
    PETScKrylovSolver(KrylovMethod method, Preconditioner pc);

    /// Create Krylov solver for a particular method and PETScPreconditioner
    PETScKrylovSolver(KrylovMethod method, PETScPreconditioner& PETScPreconditioner);

    /// Destructor
    ~PETScKrylovSolver();

    /// Solve linear system Ax = b and return number of iterations
    uint solve(const PETScMatrix& A, PETScVector& x, const PETScVector& b);
          
    /// Solve linear system Ax = b and return number of iterations
    uint solve(const PETScKrylovMatrix& A, PETScVector& x, const PETScVector& b);
    
    /// Display solver data
    void disp() const;
     
  private:

    /// Initialize KSP solver
    void init(uint M, uint N);

    /// Read parameters from database
    void readParameters();
    
    /// Set solver
    void setSolver();

    /// Set PETScPreconditioner
    void setPETScPreconditioner();
    
    /// Report the number of iterations
    void writeReport(int num_iterations);

    /// Get PETSc method identifier 
    KSPType getType(KrylovMethod method) const;

    /// Krylov method
    KrylovMethod method;

    /// PETSc PETScPreconditioner
    Preconditioner pc_petsc;

    /// DOLFIN PETScPreconditioner
    PETScPreconditioner* pc_dolfin;

    /// PETSc solver pointer
    KSP ksp;

    /// Size of old system (need to reinitialize when changing)
    uint M;
    uint N;

    /// True if we have read parameters
    bool parameters_read;
    
  };

}

#endif

#endif
