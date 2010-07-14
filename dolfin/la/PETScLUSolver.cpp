// Copyright (C) 2005-2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2009-2010.
// Modified by Niclas Jansson, 2009.
//
// First added:  2005
// Last changed: 2010-04-05

#ifdef HAS_PETSC

#include <boost/assign/list_of.hpp>
#include <dolfin/common/constants.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/main/MPI.h>
#include "LUSolver.h"
#include "PETScMatrix.h"
#include "PETScVector.h"
#include "PETScLUSolver.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
namespace dolfin
{
  class PETScKSPDeleter
  {
  public:
    void operator() (KSP* ksp)
    {
      if (ksp)
        KSPDestroy(*ksp);
      delete ksp;
    }
  };
}
//-----------------------------------------------------------------------------
// Available LU solver
const std::map<std::string, const MatSolverPackage> PETScLUSolver::lu_packages
  = boost::assign::map_list_of("default", "")
                              ("umfpack",      MAT_SOLVER_UMFPACK)
                              ("mumps",        MAT_SOLVER_MUMPS)
                              ("spooles",      MAT_SOLVER_SPOOLES)
                              ("superlu_dist", MAT_SOLVER_SUPERLU_DIST)
                              ("superlu",      MAT_SOLVER_SUPERLU);
//-----------------------------------------------------------------------------
Parameters PETScLUSolver::default_parameters()
{
  Parameters p(LUSolver::default_parameters());
  p.rename("petsc_lu_solver");
  return p;
}
//-----------------------------------------------------------------------------
PETScLUSolver::PETScLUSolver(std::string lu_package) : lu_package(lu_package)
{
  select_solver();

  // Set parameter values
  parameters = default_parameters();

  // Initialize PETSc LU solver
  init_solver();
}
//-----------------------------------------------------------------------------
PETScLUSolver::PETScLUSolver(const GenericMatrix& A, std::string lu_package)
               : lu_package(lu_package),
                 A(reference_to_no_delete_pointer(A.down_cast<PETScMatrix>()))
{
  // Check dimensions
  if (A.size(0) != A.size(1))
    error("Cannot LU factorize non-square PETSc matrix.");

  // Set parameter values
  parameters = default_parameters();

  // Select LU solver
  select_solver();

  // Initialize PETSc LU solver
  init_solver();
}
//-----------------------------------------------------------------------------
PETScLUSolver::PETScLUSolver(boost::shared_ptr<const GenericMatrix> A,
                             std::string lu_package) : lu_package(lu_package),
                             A(A)
{
  // Check dimensions
  if (A->size(0) != A->size(1))
    error("Cannot LU factorize non-square PETSc matrix.");

  // Set parameter values
  parameters = default_parameters();

  // Select LU solver
  select_solver();

  // Initialize PETSc LU solver
  init_solver();
}
//-----------------------------------------------------------------------------
PETScLUSolver::~PETScLUSolver()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void PETScLUSolver::set_operator(const GenericMatrix& A)
{
  this->A = reference_to_no_delete_pointer(A);
  assert(this->A);

  // Downcast matrix
  const PETScMatrix& _A = A.down_cast<PETScMatrix>();

  // Get some parameters
  const bool reuse_fact   = parameters["reuse_factorization"];
  const bool same_pattern = parameters["same_nonzero_pattern"];

  // Set operators with appropriate option
  if (reuse_fact)
    KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), SAME_PRECONDITIONER);
  else if (same_pattern)
    KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), SAME_NONZERO_PATTERN);
  else
    KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), DIFFERENT_NONZERO_PATTERN);

  // Check dimensions
  if (A.size(0) != A.size(1))
    error("Cannot LU factorize non-square PETSc matrix.");
}
//-----------------------------------------------------------------------------
dolfin::uint PETScLUSolver::solve(GenericVector& x, const GenericVector& b)
{
  assert(A);

  // Downcast matrix and vectors
  const PETScVector& _b = b.down_cast<PETScVector>();
  PETScVector& _x = x.down_cast<PETScVector>();

  // Check dimensions
  if (A->size(0) != b.size())
    error("Cannot LU factorize non-square PETSc matrix.");

  // Initialize solution vector (remains untouched if dimensions match)
  x.resize(A->size(1));

  // Set PETSc operators (depends on factorization re-use options);
  set_petsc_operators();

  // Write a pre-solve message
  pre_report(A->down_cast<PETScMatrix>());

  // Solve linear system
  KSPSolve(*ksp, *_b.vec(), *_x.vec());

  return 1;
}
//-----------------------------------------------------------------------------
dolfin::uint PETScLUSolver::solve(const GenericMatrix& A, GenericVector& x,
                                  const GenericVector& b)
{
  return solve(A.down_cast<PETScMatrix>(), x.down_cast<PETScVector>(),
               b.down_cast<PETScVector>());
}
//-----------------------------------------------------------------------------
dolfin::uint PETScLUSolver::solve(const PETScMatrix& A, PETScVector& x,
                                  const PETScVector& b)
{
  // Set operator
  set_operator(A);

  // Initialize solution vector (remains untouched if dimensions match)
  x.resize(A.size(1));

  // Write a pre-solve message
  pre_report(A);

  // Solve linear system
  KSPSolve(*ksp, *b.vec(), *x.vec());

  return 1;
}
//-----------------------------------------------------------------------------
std::string PETScLUSolver::str(bool verbose) const
{
  std::stringstream s;

  if (verbose)
  {
    warning("Verbose output for PETScLUSolver not implemented, calling PETSc KSPView directly.");
    KSPView(*ksp, PETSC_VIEWER_STDOUT_WORLD);
  }
  else
    s << "<PETScLUSolver>";

  return s.str();
}
//-----------------------------------------------------------------------------
void PETScLUSolver::select_solver()
{
  // Check package string
  if (lu_packages.count(lu_package) == 0)
    error("Requested PETSc LU solver '%s' is unknown,", lu_package.c_str());

  // Choose appropriate 'default' solver
  if (lu_package == "default")
  {
    if (MPI::num_processes() == 1)
      this->lu_package = "umfpack";
    else
    {
      #if PETSC_HAVE_MUMPS
      this->lu_package = "mumps";
      #elif PETSC_HAVE_SPOOLES
      this->lu_package = "spooles";
      #elif PETSC_HAVE_SUPERLU_DIST
      this->lu_package = "superlu_dist";
      #else
      error("No suitable solver for parallel LU. Consider configuring PETSc with MUMPS or SPOOLES.");
      #endif
    }
  }
}
//-----------------------------------------------------------------------------
void PETScLUSolver::init_solver()
{

  // Destroy old solver environment if necessary
  if (ksp)
  {
    if (!ksp.unique())
      error("Cannot create new KSP Krylov solver. More than one object points to the underlying PETSc object.");
  }
  ksp.reset(new KSP, PETScKSPDeleter());

  // Create solver
  if (MPI::num_processes() > 1)
    KSPCreate(PETSC_COMM_WORLD, ksp.get());
  else
    KSPCreate(PETSC_COMM_SELF, ksp.get());

  // Make solver preconditioner only
  KSPSetType(*ksp, KSPPREONLY);

  // Set preconditioner to LU factorization
  PC pc;
  KSPGetPC(*ksp, &pc);
  PCSetType(pc, PCLU);

  // Set solver package
  PCFactorSetMatSolverPackage(pc, lu_packages.find(lu_package)->second);

  // Allow matrices with zero diagonals to be solved
  #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1
  PCFactorSetShiftType(pc, MAT_SHIFT_NONZERO);
  PCFactorSetShiftAmount(pc, PETSC_DECIDE);
  #else
  PCFactorSetShiftNonzero(pc, PETSC_DECIDE);
  #endif
}
//-----------------------------------------------------------------------------
void PETScLUSolver::set_petsc_operators()
{
  // Downcast matrix
  const PETScMatrix& _A = this->A->down_cast<PETScMatrix>();

  // Get some parameters
  const bool reuse_fact   = parameters["reuse_factorization"];
  const bool same_pattern = parameters["same_nonzero_pattern"];

  // Check if operators have been set
  PetscTruth matset, pmatset;
  KSPGetOperatorsSet(*ksp, &matset, &pmatset);
  if (matset == PETSC_TRUE && pmatset == PETSC_TRUE)
  {
    // Get preconditioner re-use flag
    MatStructure mat_struct;
    Mat Amat, Pmat;
    KSPGetOperators(*ksp, &Amat, &Pmat, &mat_struct);

    // Set operators if necessary
    if (reuse_fact && mat_struct != SAME_PRECONDITIONER)
      KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), SAME_PRECONDITIONER);
    else if (same_pattern && mat_struct != SAME_NONZERO_PATTERN)
      KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), SAME_NONZERO_PATTERN);
    else if (!reuse_fact && !same_pattern &&  mat_struct != DIFFERENT_NONZERO_PATTERN)
      KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), DIFFERENT_NONZERO_PATTERN);
  }
  else if (matset != PETSC_TRUE && pmatset != PETSC_TRUE)
  {
    // Set operators with appropriate option
    if (reuse_fact)
      KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), SAME_PRECONDITIONER);
    else if (same_pattern)
      KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), SAME_NONZERO_PATTERN);
    else
      KSPSetOperators(*ksp, *_A.mat(), *_A.mat(), DIFFERENT_NONZERO_PATTERN);
  }
  else
    error("Operators incorrectly set for PETSc.");
}
//-----------------------------------------------------------------------------
void PETScLUSolver::pre_report(const PETScMatrix& A) const
{
  const MatSolverPackage solver_type;
  PC pc;
  KSPGetPC(*ksp, &pc);
  PCFactorGetMatSolverPackage(pc, &solver_type);

  // Get parameter
  const bool report = parameters["report"];

  if (report)
  {
    info(PROGRESS, "Solving linear system of size %d x %d (PETSc LU solver, %s).",
         A.size(0), A.size(1), solver_type);
  }
}
//-----------------------------------------------------------------------------

#endif
