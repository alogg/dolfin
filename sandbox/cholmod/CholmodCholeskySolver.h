// Copyright (C) 2008 Dag Lindbo and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-08-15
// Last changed: 2008-08-15

#ifndef __CHOLMOD_CHOLESKY_SOLVER_H
#define __CHOLMOD_CHOLESKY_SOLVER_H

#include <dolfin/parameter/Parametrized.h>

extern "C" 
{
#ifdef HAS_UMFPACK
#include <cholmod.h>
#endif
}

namespace dolfin
{
  /// Forward declarations
  class GenericVector;
  class GenericMatrix;

  /// This class implements the direct solution (Cholesky factorization) of
  /// linear systems of the form Ax = b. Sparse matrices
  /// are solved using CHOLMOD http://www.cise.ufl.edu/research/sparse/cholmod/
  /// is installed.
    
  class CholmodCholeskySolver : public Parametrized
  {

  public:
    
    /// Constructor
    CholmodCholeskySolver();

    /// Destructor
    ~CholmodCholeskySolver();

    /// Solve linear system Ax = b for a sparse matrix using CHOLMOD
    virtual uint solve(const GenericMatrix& A, GenericVector& x, const GenericVector& b);

    /// Cholesky-factor sparse matrix A if CHOLMOD is installed
    virtual uint factorize(const GenericMatrix& A);

    /// Solve factorized system (CHOLMOD).
    virtual uint factorizedSolve(GenericVector& x, const GenericVector& b);

  private:
    
    // Data for Cholesky factorization of sparse ublas matrix (cholmod only)
    class Cholmod
    {
    public:
 
      Cholmod();
      ~Cholmod();

      // Clear data
      void clear();

      // Initialise with matrix
      void init(long int* Ap, long int* Ai, double* Ax, uint M, uint nz);

      // Factorize
      void factorize();

      // Factorized solve
      void factorizedSolve(double*x, const double* b);

      /// Check status flag returned by an CHOLMOD function
      void checkStatus(long int status, std::string function) const;

      // CHOLMOD data
      cholmod_sparse *A_chol;
      cholmod_factor *L_chol;
      cholmod_common c;

      uint N;
      bool factorized;
    };

    Cholmod cholmod;
  };

}

#endif
