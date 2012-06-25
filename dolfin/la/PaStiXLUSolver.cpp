// Copyright (C) 2011 Garth N. Wells
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2011-10-19
// Last changed:

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include "dolfin/common/MPI.h"
#include "dolfin/common/NoDeleter.h"
#include "dolfin/common/types.h"
#include "dolfin/log/dolfin_log.h"
#include "dolfin/parameter/GlobalParameters.h"
#include "GenericVector.h"
#include "SparsityPattern.h"
#include "GenericMatrix.h"
#include "LUSolver.h"
#include "STLMatrix.h"
#include "PaStiXLUSolver.h"

#ifdef HAS_PASTIX

extern "C"
{
#include <pastix.h>
}

using namespace dolfin;

//-----------------------------------------------------------------------------
Parameters PaStiXLUSolver::default_parameters()
{
  Parameters p(LUSolver::default_parameters());
  p.rename("pastix_lu_solver");

  // Number of threads per MPI process
  p.add<uint>("num_threads");

  // Check matrix for consistency
  p.add("check_matrix", false);

  return p;
}
//-----------------------------------------------------------------------------
PaStiXLUSolver::PaStiXLUSolver(const STLMatrix& A)
  : A(reference_to_no_delete_pointer(A)), id(0)
{
  // Set parameter values
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
PaStiXLUSolver::PaStiXLUSolver(boost::shared_ptr<const STLMatrix> A)
  : A(A), id(0)
{
  // Set parameter values
  parameters = default_parameters();
}
//-----------------------------------------------------------------------------
unsigned int PaStiXLUSolver::solve(GenericVector& x, const GenericVector& b)
{
  assert(A);

  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  int    iparm[IPARM_SIZE];
  double dparm[DPARM_SIZE];
  for (int i = 0; i < IPARM_SIZE; i++)
    iparm[i] = 0;
  for (int i = 0; i < DPARM_SIZE; i++)
    dparm[i] = 0;

  // Set default parameters
  pastix_initParam(iparm, dparm);

  // LU or Cholesky
  const bool symmetric = parameters["symmetric_operator"];
  if (symmetric)
  {
    iparm[IPARM_SYM] = API_SYM_YES;
    iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
  }
  else
  {
    iparm[IPARM_SYM] = API_SYM_NO;
    iparm[IPARM_FACTORIZATION] = API_FACT_LU;
  }

  // Matrix data in compressed sparse column format (C indexing)
  std::vector<double> vals;
  std::vector<uint> rows, col_ptr, local_to_global_cols;
  A->csc(vals, rows, col_ptr, local_to_global_cols, symmetric);

  dolfin_assert(local_to_global_cols.size() > 0);

  const std::vector<uint> local_to_global_cols_ref = local_to_global_cols;

  int* _col_ptr = reinterpret_cast<int*>(col_ptr.data());
  int* _rows = reinterpret_cast<int*>(rows.data());
  int* _local_to_global_cols = reinterpret_cast<int*>(local_to_global_cols.data());
  double* _vals = vals.data();

  const uint n = col_ptr.size() - 1;

  // Check matrix
  if (parameters["check_matrix"])
  {
    d_pastix_checkMatrix(mpi_comm, API_VERBOSE_YES, iparm[IPARM_SYM], API_YES,
  		                   n, &_col_ptr, &_rows, &_vals, &_local_to_global_cols, 1);
  }

  // PaStiX object
  pastix_data_t* pastix_data = NULL;

  // Number of threads per MPI process
  if (parameters["num_threads"].is_set())
    iparm[IPARM_THREAD_NBR] = std::max((uint) 1, (uint) parameters["num_threads"]);
  else
    iparm[IPARM_THREAD_NBR] = std::max((uint) 1, (uint) dolfin::parameters["num_threads"]);

  // User-supplied RHS
  iparm[IPARM_RHS_MAKING] = API_RHS_B;

  // Level of verbosity
  if (parameters["verbose"])
    iparm[IPARM_VERBOSE] = API_VERBOSE_YES;
  else
    iparm[IPARM_VERBOSE] = API_VERBOSE_NO;

  // Graph (matrix) distributed
  if (MPI::num_processes() > 1)
    iparm[IPARM_GRAPHDIST] = API_YES;
  else
    iparm[IPARM_GRAPHDIST] = API_NO;

  dolfin_assert(local_to_global_cols.size() > 0);
  std::vector<int> perm(local_to_global_cols.size());
  std::vector<int> invp(local_to_global_cols.size());

  // Number of RHS vectors
  const int nrhs = 1;

  // Re-order
  iparm[IPARM_START_TASK] = API_TASK_ORDERING;
  iparm[IPARM_END_TASK]   = API_TASK_BLEND;
  d_dpastix(&pastix_data, mpi_comm, n, _col_ptr, _rows, _vals,
            _local_to_global_cols,
            perm.data(), invp.data(),
            NULL, nrhs, iparm, dparm);

  // Factorise
  iparm[IPARM_START_TASK] = API_TASK_NUMFACT;
  iparm[IPARM_END_TASK]   = API_TASK_NUMFACT;
  d_dpastix(&pastix_data, mpi_comm, n, _col_ptr, _rows, _vals,
            _local_to_global_cols,
            perm.data(), invp.data(),
            NULL, nrhs, iparm, dparm);

  // Get RHS data for this process
  std::vector<double> _b;
  b.gather(_b, local_to_global_cols_ref);
  double* b_ptr = _b.data();

  // Solve
  iparm[IPARM_START_TASK] = API_TASK_SOLVE;
  iparm[IPARM_END_TASK] = API_TASK_SOLVE;
  d_dpastix(&pastix_data, mpi_comm, n, _col_ptr, _rows, _vals,
            _local_to_global_cols,
            perm.data(), invp.data(),
            b_ptr, nrhs, iparm, dparm);

  // Distribute solution
  assert(b.size() == x.size());
  x.set(_b.data(), local_to_global_cols_ref.size(), local_to_global_cols_ref.data());
  x.apply("insert");

  // Clean up
  iparm[IPARM_START_TASK] = API_TASK_CLEAN;
  iparm[IPARM_END_TASK] = API_TASK_CLEAN;
  d_dpastix(&pastix_data, mpi_comm, n, NULL, NULL, NULL,
            _local_to_global_cols,
            perm.data(), invp.data(),
            b_ptr, nrhs, iparm, dparm);

  return 1;
}
//-----------------------------------------------------------------------------
PaStiXLUSolver::~PaStiXLUSolver()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
#endif
