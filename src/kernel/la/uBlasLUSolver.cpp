// Copyright (C) 2006-2007 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg 2006.
// 
// First added:  2006-06-01
// Last changed: 2007-07-13

#include <dolfin/dolfin_log.h>
#include <dolfin/uBlasLUSolver.h>
#include <dolfin/uBlasKrylovSolver.h>
#include <dolfin/uBlasSparseMatrix.h>
#include <dolfin/uBlasKrylovMatrix.h>
#include <dolfin/uBlasVector.h>

extern "C" 
{
// Take care of different default locations
#ifdef HAVE_UMFPACK_H
  #include <umfpack.h>
#elif HAVE_UMFPACK_UMFPACK_H
  #include <umfpack/umfpack.h>
#elif HAVE_UFSPARSE_UMFPACK_H
  #include <ufsparse/umfpack.h>
#elif HAVE_SUITESPARSE_UMFPACK_H
  #include <suitesparse/umfpack.h>
#endif
}

using namespace dolfin;

//-----------------------------------------------------------------------------
uBlasLUSolver::uBlasLUSolver() : uBlasLinearSolver(), AA(0), ej(0), Aj(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
uBlasLUSolver::~uBlasLUSolver()
{
  if ( AA ) delete AA;
  if ( ej ) delete ej;
  if ( Aj ) delete Aj;
}
//-----------------------------------------------------------------------------
dolfin::uint uBlasLUSolver::solve(const uBlasMatrix<ublas_dense_matrix>& A, 
                                  uBlasVector& x, const uBlasVector& b)
{    
  // Make copy of matrix and vector
  ublas_dense_matrix Atemp(A);
  x.resize(b.size());
  x.assign(b);

  // Solve
  return solveInPlace(Atemp, x);
}
//-----------------------------------------------------------------------------
#if defined(HAVE_UMFPACK_H) || defined(HAVE_UMFPACK_UMFPACK_H) || defined(HAVE_UFSPARSE_UMFPACK_H) ||  defined(HAVE_SUITESPARSE_UMFPACK_H)
dolfin::uint uBlasLUSolver::solve(const uBlasMatrix<ublas_sparse_matrix>& A, uBlasVector& x, 
                                  const uBlasVector& b)
{
  // Check dimensions and get number of non-zeroes
  const uint M  = A.size(0);
  const uint N  = A.size(1);
  const uint nz = A.nnz();

  dolfin_assert(M == A.size(1));
  dolfin_assert(nz > 0);

  x.init(N);

  // Make sure matrix assembly is complete
  (const_cast< uBlasMatrix<ublas_sparse_matrix>& >(A)).complete_index1_data(); 

  message("Solving linear system of size %d x %d (UMFPACK LU solver).", M, N);

  //FIXME: From UMFPACK v.5.0 onwards, UF_long is introduced and should be used 
  //       in place of long int.

  double* dnull = (double *) NULL;
  long int*    inull = (long int *) NULL;
  void *Symbolic, *Numeric;
  const std::size_t* Ap = &(A.index1_data() [0]);
  const std::size_t* Ai = &(A.index2_data() [0]);
  const double* Ax = &(A.value_data() [0]);
  double* xx = &(x.data() [0]);
  const double* bb = &(b.data() [0]);

  // Solve for transpose since we use compressed row format, and UMFPACK 
  // expects compressed column format

  long int* Rp = new long int[M+1];
  long int* Ri = new long int[nz];
  double* Rx   = new double[nz];

  long int status;

  // Compute transpose
  status= umfpack_dl_transpose(M, M, (const long int*) Ap, (const long int*) Ai, Ax, inull, inull, Rp, Ri, Rx);
  check_status(status, "transpose");

  // Solve procedure
  status= umfpack_dl_symbolic(M, M, (const long int*) Rp, (const long int*) Ri, Rx, &Symbolic, dnull, dnull);
  check_status(status, "symbolic");

  status = umfpack_dl_numeric( (const long int*) Rp, (const long int*) Ri, Rx, Symbolic, &Numeric, dnull, dnull);
  check_status(status, "numeric");

  umfpack_dl_free_symbolic(&Symbolic);

  status = umfpack_dl_solve(UMFPACK_A, (const long int*) Rp, (const long int*) Ri, Rx, xx, bb, Numeric, dnull, dnull);
  check_status(status, "solve");
 
  umfpack_dl_free_numeric(&Numeric);

  // Clean up
  delete [] Rp;
  delete [] Ri;
  delete [] Rx;

  return 1;
}

#else

dolfin::uint uBlasLUSolver::solve(const uBlasMatrix<ublas_sparse_matrix>& A, uBlasVector& x, 
    const uBlasVector& b)
{
  warning("UMFPACK must be installed to peform a LU solve for uBlas matrices. A Krylov iterative solver will be used instead.");

  uBlasKrylovSolver solver;
  return solver.solve(A, x, b);
}
#endif
//-----------------------------------------------------------------------------
void uBlasLUSolver::solve(const uBlasKrylovMatrix& A, uBlasVector& x,
			  const uBlasVector& b)
{
  // The linear system is solved by computing a dense copy of the matrix,
  // obtained through multiplication with unit vectors.

  // Check dimensions
  const uint M  = A.size(0);
  const uint N  = A.size(1);
  dolfin_assert(M == N);
  dolfin_assert(M == b.size());

  // Initialize temporary data if not already done
  if ( !AA )
  {
    AA = new uBlasMatrix<ublas_dense_matrix>(M, N);
    ej = new uBlasVector(N);
    Aj = new uBlasVector(M);
  }
  else
  {
    AA->init(M, N);
    ej->init(N);
    Aj->init(N);
  }

  // Reset unit vector
  *ej = 0.0;

  // Compute columns of matrix
  for (uint j = 0; j < N; j++)
  {
    (*ej)(j) = 1.0;

    // Compute product Aj = Aej
    A.mult(*ej, *Aj);
    
    // Set column of A
    column(*AA, j) = *Aj;
    
    (*ej)(j) = 0.0;
  }

  // Solve linear system
  solve(*AA, x, b);
}
//-----------------------------------------------------------------------------
dolfin::uint uBlasLUSolver::solveInPlaceUBlas(uBlasMatrix<ublas_dense_matrix>& A, 
                                      uBlasVector& x, const uBlasVector& b) const
{
  const uint M = A.size1();
  dolfin_assert(M == b.size());
  
  if( x.size() != M )
    x.resize(M);

  // Initialise solution vector
  x.assign(b);

  // Solve
  return solveInPlace(A, x);
}
//-----------------------------------------------------------------------------
void uBlasLUSolver::check_status(long int status, std::string function) const
{
#if defined(HAVE_UMFPACK_H)|| defined(HAVE_UMFPACK_UMFPACK_H) || defined(HAVE_UFSPARSE_UMFPACK_H) ||  defined(HAVE_SUITESPARSE_UMFPACK_H)
  if(status == UMFPACK_OK)
    return;

  // Help out by printing which UMFPACK function is returning the warning/error
  cout << "UMFPACK problem related to call to " << function << endl;

  if(status == UMFPACK_WARNING_singular_matrix)
    warning("UMFPACK reports that the matrix being solved is singular.");
  else if(status == UMFPACK_ERROR_out_of_memory)
    error("UMFPACK has run out of memory solving a system.");
  else if(status == UMFPACK_ERROR_invalid_system)
    error("UMFPACK reports an invalid system. Is the matrix square?.");
  else if(status == UMFPACK_ERROR_invalid_Numeric_object)
    error("UMFPACK reports an invalid Numeric object.");
  else if(status == UMFPACK_ERROR_invalid_Symbolic_object)
    error("UMFPACK reports an invalid Symbolic object.");
  else if(status != UMFPACK_OK)
    warning("UMFPACK is reporting an unknown error.");

#else

  error("Problem with DOLFIN build configuration for using UMFPACK.");   

#endif
}
//-----------------------------------------------------------------------------



