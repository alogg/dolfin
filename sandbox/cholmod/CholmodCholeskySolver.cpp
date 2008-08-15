// Copyright (C) 2006-2008 Dag Lindbo and Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-08-15
// Last changed: 2008-08-15

#include <dolfin/log/dolfin_log.h>
#include "CholmodCholeskySolver.h"
#include <dolfin/la/GenericMatrix.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/KrylovSolver.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
CholmodCholeskySolver::CholmodCholeskySolver()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
CholmodCholeskySolver::~CholmodCholeskySolver()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
#ifdef HAS_UMFPACK
dolfin::uint CholmodCholeskySolver::solve(const GenericMatrix& A,
					  GenericVector& x, 
					  const GenericVector& b)
{
  // Factorize matrix
  factorize(A);

  // Solve system
  factorizedSolve(x, b);

  // Clear data
  cholmod.clear();

  return 1;
}
//-----------------------------------------------------------------------------
dolfin::uint CholmodCholeskySolver::factorize(const GenericMatrix& A)
{
  // Check dimensions and get number of non-zeroes
  boost::tuple<const std::size_t*, const std::size_t*, const double*, int> data = A.data();
  const uint M   = A.size(0);
  const uint nnz = data.get<3>();
  dolfin_assert(A.size(0) == A.size(1));

  dolfin_assert(nnz >= M); 

  // Initialise cholmod data
  // NOTE: Casting away const here
  cholmod.init((UF_long*) data.get<0>(),(UF_long*) data.get<1>(), 
	       (double*) data.get<2>(), M, nnz);

  // Factorize
  message("Cholesky-factorizing linear system of size %d x %d (CHOLMOD).",M,M);
  cholmod.factorize();

  return 1;
}
//-----------------------------------------------------------------------------
dolfin::uint CholmodCholeskySolver::factorizedSolve(GenericVector& x, const GenericVector& b)
{
  const uint N = b.size();

  if(!cholmod.factorized)
    error("Factorized solve must be preceeded by call to factorize.");

  if(N != cholmod.N)
    error("Vector does not match size of factored matrix");

  // Initialise solution vector and solve
  x.init(N);

  message("Solving factorized linear system of size %d x %d (CHOLMOD).", N, N);

  cholmod.factorizedSolve(x.data(), b.data());

  return 1;
}
//-----------------------------------------------------------------------------
#else
dolfin::uint CholmodCholeskySolver::solve(const GenericMatrix& A, 
					  GenericVector& x, 
					  const GenericVector& b)
{
  warning("CHOLMOD must be installed to peform a Cholesky solve for the current backend. A Krylov iterative solver will be used instead.");

  KrylovSolver solver;
  return solver.solve(A, x, b);
}
//-----------------------------------------------------------------------------
dolfin::uint CholmodCholeskySolver::factorize(const GenericMatrix& A)
{
  error("CHOLMOD must be installed to perform sparse Cholesky factorization.");
  return 0;
}
//-----------------------------------------------------------------------------
dolfin::uint CholmodCholeskySolver::factorizedSolve(GenericVector& x, const GenericVector& b) const
{
  error("CHOLMOD must be installed to perform sparse back and forward substitution");
  return 0;
}
#endif
//-----------------------------------------------------------------------------
// CholmodCholeskySolver::Cholmod implementation
//-----------------------------------------------------------------------------
CholmodCholeskySolver::Cholmod::Cholmod() : 
  A_chol(0), L_chol(0), N(0), factorized(false) 
{ 
  // "Start" cholmod
  cholmod_l_start(&c);
}
//-----------------------------------------------------------------------------
CholmodCholeskySolver::Cholmod::~Cholmod() 
{ 
  clear(); 
	
  // "stop" cholmod
  cholmod_l_finish(&c);
}
//-----------------------------------------------------------------------------
void CholmodCholeskySolver::Cholmod::clear()
{

  if(A_chol)
    {
      cholmod_l_free(1, sizeof(cholmod_sparse), A_chol, &c);
      A_chol = 0;
    }
  if(L_chol)
    {
      cholmod_l_free_factor(&L_chol, &c);
      L_chol = 0;
    }
}
//-----------------------------------------------------------------------------
void CholmodCholeskySolver::Cholmod::init(UF_long* Ap, 
					  UF_long* Ai, 
					  double* Ax, 
					  uint M, uint nz)
{  
  if(factorized)
    warning("CholeskySolver already contains a factorized matrix! Clearing and starting over.");

  // Clear any data
  clear();

  A_chol = (cholmod_sparse*) cholmod_l_malloc(1, sizeof(cholmod_sparse), &c);

  // The matrix
  A_chol->p = Ap;
  A_chol->i = Ai;
  A_chol->x = Ax;

  A_chol->nzmax = nz;
  A_chol->ncol = M;
  A_chol->nrow = M;

  A_chol->sorted = 1;
  A_chol->packed = 1;

  A_chol->xtype = CHOLMOD_REAL;
  A_chol->stype = -1;
  A_chol->dtype = CHOLMOD_DOUBLE;
  A_chol->itype = CHOLMOD_LONG;

  N = M;
}
//-----------------------------------------------------------------------------
void CholmodCholeskySolver::Cholmod::factorize()
{
#ifdef HAS_UMFPACK

  // Analyze
  L_chol = cholmod_l_analyze(A_chol, &c);

  // Factorize
  cholmod_l_factorize(A_chol, L_chol, &c);

  factorized = true;

#else
  error("CHOLMOD not installed");
#endif
}
//-----------------------------------------------------------------------------
void CholmodCholeskySolver::Cholmod::factorizedSolve(double*x, const double* b)
{
#ifdef HAS_UMFPACK

  cholmod_dense *x_chol, *b_chol;

  // initialize rhs
  b_chol = (cholmod_dense*) cholmod_l_malloc(1, sizeof(cholmod_dense), &c);
  b_chol->x = (double*) b;
  b_chol->nrow = N;
  b_chol->ncol = 1;
  b_chol->d = N;
  b_chol->nzmax = N;
  b_chol->xtype = A_chol->xtype;
  b_chol->dtype = A_chol->dtype;
  
  // solve
  x_chol = cholmod_l_solve(CHOLMOD_A, L_chol, b_chol, &c);

  // solution vector
  // FIXME
  memcpy(x, x_chol->x, N*sizeof(double));
  cholmod_l_free_dense(&x_chol, &c);

  // clear rhs
  cholmod_l_free(1, sizeof(cholmod_dense), b_chol, &c);

#else
  error("CHOLMOD not installed");
#endif
}
//-----------------------------------------------------------------------------
void CholmodCholeskySolver::Cholmod::checkStatus(UF_long status, std::string function) const
{
#ifdef HAS_UMFPACK
  if(status == CHOLMOD_OK)
    return;

  // Printing which CHOLMOD function is returning an warning/error
  cout << "CHOLMOD problem related to call to " << function << endl;

  // if(status == CHOLMOD_WARNING_singular_matrix)
  //   warning("CHOLMOD reports that the matrix being solved is singular.");
  // else if(status == CHOLMOD_ERROR_out_of_memory)
  //   error("CHOLMOD has run out of memory solving a system.");
  // else if(status == CHOLMOD_ERROR_invalid_system)
  //   error("CHOLMOD reports an invalid system. Is the matrix square?.");
  // else if(status == CHOLMOD_ERROR_invalid_Numeric_object)
  //   error("CHOLMOD reports an invalid Numeric object.");
  // else if(status == CHOLMOD_ERROR_invalid_Symbolic_object)
  //   error("CHOLMOD reports an invalid Symbolic object.");
  // else if(status != CHOLMOD_OK)
  //   warning("CHOLMOD is reporting an unknown error.");

#else
  error("Problem with DOLFIN build configuration for using CHOLMOD.");   
#endif
}
//-----------------------------------------------------------------------------
