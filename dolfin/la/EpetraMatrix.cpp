// Copyright (C) 2008 Martin Sandve Alnes, Kent-Andre Mardal and Johannes Ring.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
// Modified by Garth N. Wells, 2008.
//
// First added:  2008-04-21
// Last changed: 2008-12-30

#ifdef HAS_TRILINOS

#include <cstring>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/common/Array.h>
#include "EpetraVector.h"
#include "EpetraMatrix.h"
#include "GenericSparsityPattern.h"
#include "EpetraSparsityPattern.h"
#include "EpetraFactory.h"
#include <dolfin/common/Timer.h>

#include <Epetra_CrsGraph.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <EpetraExt_MatrixMatrix.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
EpetraMatrix::EpetraMatrix():
    Variable("A", "Epetra matrix"),
    A(static_cast<Epetra_FECrsMatrix*>(0))
{
  // TODO: call Epetra_Init or something?
}
//-----------------------------------------------------------------------------
EpetraMatrix::EpetraMatrix(uint M, uint N):
    Variable("A", "Epetra matrix"),
    A(static_cast<Epetra_FECrsMatrix*>(0))
{
  // TODO: call Epetra_Init or something?
  // Create Epetra matrix
  resize(M, N);
}
//-----------------------------------------------------------------------------
EpetraMatrix::EpetraMatrix(const EpetraMatrix& A):
  Variable("A", "Epetra matrix"),
  A(static_cast<Epetra_FECrsMatrix*>(0))
{
  if (A.mat())
  {
    std::tr1::shared_ptr<Epetra_FECrsMatrix> _A(new Epetra_FECrsMatrix(*A.mat()));
    this->A = _A;
  }
}
//-----------------------------------------------------------------------------
EpetraMatrix::EpetraMatrix(std::tr1::shared_ptr<Epetra_FECrsMatrix> A):
    Variable("A", "a sparse matrix"),
    A(A)
{
  // TODO: call Epetra_Init or something?
}
//-----------------------------------------------------------------------------
EpetraMatrix::EpetraMatrix(const Epetra_CrsGraph& graph):
    Variable("A", "a sparse matrix"),
    A(new Epetra_FECrsMatrix(Copy, graph))
{
  // TODO: call Epetra_Init or something?
}
//-----------------------------------------------------------------------------
EpetraMatrix::~EpetraMatrix()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void EpetraMatrix::resize(uint M, uint N)
{
  // Not yet implemented
  error("EpetraMatrix::resize(uint, unit) not yet implemented.");
}
//-----------------------------------------------------------------------------
void EpetraMatrix::init(const GenericSparsityPattern& sparsity_pattern)
{
  if (!A.unique())
    error("Cannot initialise EpetraMatrix. More than one object points to the underlying Epetra object.");

  const EpetraSparsityPattern& epetra_pattern = dynamic_cast<const EpetraSparsityPattern&>(sparsity_pattern);
  std::tr1::shared_ptr<Epetra_FECrsMatrix> _A(new Epetra_FECrsMatrix(Copy, epetra_pattern.pattern()));
  A = _A;
}
//-----------------------------------------------------------------------------
EpetraMatrix* EpetraMatrix::copy() const
{
  EpetraMatrix* mcopy = new EpetraMatrix(*this);
  return mcopy;
}
//-----------------------------------------------------------------------------
dolfin::uint EpetraMatrix::size(uint dim) const
{
  dolfin_assert(A); 
  int M = A->NumGlobalRows();
  int N = A->NumGlobalCols();
  return (dim == 0 ? M : N);
}
//-----------------------------------------------------------------------------
void EpetraMatrix::get(double* block,
		       uint m, const uint* rows,
		       uint n, const uint* cols) const
{
  dolfin_assert(A);
  
  int max_num_indices = A->MaxNumEntries();
  int num_entities    = 0;
  int * indices       = new int[max_num_indices];
  double * values     = new double[max_num_indices];
  
  // For each row in rows
  for(uint i = 0; i < m; ++i)
  {
    // Extract the values and indices from row: rows[i]
    if (A->IndicesAreLocal())
      A->ExtractMyRowView(rows[i], num_entities, values, indices);
    else
      A->ExtractGlobalRowView(rows[i], num_entities, values, indices);

    int k = 0;
    // Step the indices to the start of cols
    while (indices[k] < static_cast<int>(cols[0]))
      k++;
    // Fill the collumns in the block 
    for (uint j = 0; j < n; j++)
    {
      if (k < num_entities and indices[k] == static_cast<int>(cols[j]))
      {
        block[i*n + j] = values[k];
        k++;
      }
      else
        block[i*n + j] = 0.0;
    }
  }
}
//-----------------------------------------------------------------------------
void EpetraMatrix::set(const double* block,
                       uint m, const uint* rows,
                       uint n, const uint* cols)
{
  dolfin_assert(A); 
  int err = A->ReplaceGlobalValues(m, reinterpret_cast<const int*>(rows),
                                   n, reinterpret_cast<const int*>(cols), block);
  if (err!= 0) 
    error("Did not manage to set the values into the matrix"); 
}
//-----------------------------------------------------------------------------
void EpetraMatrix::add(const double* block,
                       uint m, const uint* rows,
                       uint n, const uint* cols)
{
  Timer t0("Matrix add"); 
  dolfin_assert(A); 

  int err = A->SumIntoGlobalValues(m, reinterpret_cast<const int*>(rows), 
                                   n, reinterpret_cast<const int*>(cols), block, 
                                   Epetra_FECrsMatrix::ROW_MAJOR);

  if (err!= 0) 
    error("Did not manage to put the values into the matrix"); 
}
//-----------------------------------------------------------------------------
void EpetraMatrix::axpy(double a, const GenericMatrix& A)
{
  const EpetraMatrix* AA = &A.down_cast<EpetraMatrix>();

  dolfin_assert((*AA->mat()).NumGlobalNonzeros() == this->A->NumGlobalNonzeros() 
             && (*AA->mat()).NumGlobalCols()     == this->A->NumGlobalCols() 
             && (*AA->mat()).NumGlobalRows()     == this->A->NumGlobalRows() );

  EpetraExt::MatrixMatrix::Add(*(AA->mat()), false, a, *(this->A), 1.0);
}
//-----------------------------------------------------------------------------
void EpetraMatrix::zero()
{
  dolfin_assert(A);
  A->PutScalar(0.0);
}
//-----------------------------------------------------------------------------
void EpetraMatrix::apply()
{
  dolfin_assert(A); 
  A->GlobalAssemble();
  // TODO
  //A->OptimizeStorage(); 
}
//-----------------------------------------------------------------------------
void EpetraMatrix::disp(uint precision) const
{
  dolfin_assert(A); 
  A->Print(std::cout); 
}
//-----------------------------------------------------------------------------
void EpetraMatrix::ident(uint m, const uint* rows)
{
  dolfin_assert(A); 
  double* values;
  int* indices; 
  int* row_size = new int; 
  int r;

  for (uint i=0; i<m; i++)
  {
    r = rows[i]; 
    int err = A->ExtractMyRowView(r, *row_size, values, indices); 
    if (err!= 0) 
      error("Trouble with ExtractMyRowView in EpetraMatrix::ident."); 
    memset(values, 0, (*row_size)*sizeof(double)); 
    for (uint j=0; j<m; j++) 
    {
      if (r == indices[j]) 
      {
        values[j] = 1.0; 
        break; 
      }
    }
  }
  delete row_size; 
}
//-----------------------------------------------------------------------------
void EpetraMatrix::zero(uint m, const uint* rows)
{
  dolfin_assert(A); 
  double* values; 
  int* indices; 
  int* row_size = new int; 
  int r;

  for (uint i=0; i<m; i++)
  {
    r = rows[i]; 
    int err = A->ExtractMyRowView(r, *row_size, values, indices); 
    if (err!= 0) 
      error("Trouble with ExtractMyRowView in EpetraMatrix::ident."); 
    memset(values, 0, (*row_size)*sizeof(double)); 
  }
  delete row_size; 
}
//-----------------------------------------------------------------------------
void EpetraMatrix::mult(const GenericVector& x_, GenericVector& Ax_, bool transposed) const
{
  dolfin_assert(A); 

  const EpetraVector* x = dynamic_cast<const EpetraVector*>(x_.instance());  
  if (!x) 
    error("The vector x should be of type EpetraVector");  

  EpetraVector* Ax = dynamic_cast<EpetraVector*>(Ax_.instance());  
  if (!Ax) 
    error("The vector Ax should be of type EpetraVector");  

  A->Multiply(transposed, *(x->vec()), *(Ax->vec()));
}
//-----------------------------------------------------------------------------
void EpetraMatrix::getrow(uint row, Array<uint>& columns, Array<double>& values) const
{
  dolfin_assert(A); 

  // Temporary variables
  int *indices; 
  double* vals; 
  int* num_entries = new int; 

  // Extract data from Epetra matrix 
  int err = A->ExtractMyRowView(row, *num_entries, vals, indices); 
  if (err!= 0) 
    error("Did not manage to get a copy of the row."); 

  // Put data in columns and values
  columns.clear();
  values.clear(); 
  for (int i=0; i< *num_entries; i++)
  {
    columns.push_back(indices[i]);
    values.push_back(vals[i]);
  }

  delete num_entries; 
}
//-----------------------------------------------------------------------------
void EpetraMatrix::setrow(uint row, const Array<uint>& columns, const Array<double>& values)
{
  error("Not implemented.");
}
//-----------------------------------------------------------------------------
LinearAlgebraFactory& EpetraMatrix::factory() const
{
  return EpetraFactory::instance();
}
//-----------------------------------------------------------------------------
Epetra_FECrsMatrix& EpetraMatrix::mat_ref() const
{
  dolfin_assert(A); 
  warning("This may be removed shortly. If possible, use EpetraVector::vec() instead.");
  return *A;
}
//-----------------------------------------------------------------------------
std::tr1::shared_ptr<Epetra_FECrsMatrix> EpetraMatrix::mat() const
{
  dolfin_assert(A); 
  return A;
}
//-----------------------------------------------------------------------------
const EpetraMatrix& EpetraMatrix::operator*= (double a)
{
  dolfin_assert(A);
  A->Scale(a);
  return *this;
}
//-----------------------------------------------------------------------------
const EpetraMatrix& EpetraMatrix::operator/= (double a)
{
  dolfin_assert(A);
  A->Scale(1.0/a);
  return *this;
}
//-----------------------------------------------------------------------------
const GenericMatrix& EpetraMatrix::operator= (const GenericMatrix& A)
{
  *this = A.down_cast<EpetraMatrix>();
  return *this;
}
//-----------------------------------------------------------------------------
const EpetraMatrix& EpetraMatrix::operator= (const EpetraMatrix& A)
{
  dolfin_assert(A.mat());
  *(this->A) = *A.mat();
  return *this;
}
//-----------------------------------------------------------------------------
LogStream& dolfin::operator<< (LogStream& stream, const EpetraMatrix& A)
{
  stream << "[ Epetra matrix of size " << A.size(0) << " x " << A.size(1) << " ]";
  return stream;
}
//-----------------------------------------------------------------------------
#endif
