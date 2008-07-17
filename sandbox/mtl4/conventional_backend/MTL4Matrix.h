// Copyright (C) 2008 Dag Lindbo
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-07-06
// Last changed: 2008-07-06

#ifndef __MTL4_MATRIX_H
#define __MTL4_MATRIX_H

#ifdef HAS_MTL4

#include <dolfin/log/LogStream.h>
#include <dolfin/common/Variable.h>
#include "mtl4.h"
#include "GenericMatrix.h"

namespace dolfin
{

  class MTL4Matrix: public GenericMatrix, public Variable
  {
  public:

    /// Create empty matrix
    MTL4Matrix();

    /// Create M x N matrix
    MTL4Matrix(uint M, uint N);

    /// Copy constuctor
    explicit MTL4Matrix(const MTL4Matrix& A);

    /// Create view from given MTL4_sparse_matrix pointer
    explicit MTL4Matrix(MTL4_sparse_matrix* A);

    /// Destructor
    virtual ~MTL4Matrix();

    //--- Implementation of the GenericTensor interface ---

    /// Initialize zero tensor using sparsity pattern
    virtual void init(const GenericSparsityPattern& sparsity_pattern);

    /// Return copy of tensor
    virtual MTL4Matrix* copy() const;

    /// Return size of given dimension
    virtual uint size(uint dim) const;

    /// Set all entries to zero and keep any sparse structure
    virtual void zero();

    /// Finalize assembly of tensor
    virtual void apply();

    /// Display tensor
    virtual void disp(uint precision=2) const;

    //--- Implementation of the GenericMatrix interface ---

    /// Initialize M x N matrix
    virtual void init(uint M, uint N);

    /// Get block of values
    virtual void get(real* block, uint m, const uint* rows, uint n, const uint* cols) const;

    /// Set block of values
    virtual void set(const real* block, uint m, const uint* rows, uint n, const uint* cols);

    /// Add block of values
    virtual void add(const real* block, uint m, const uint* rows, uint n, const uint* cols);

    /// Get non-zero values of given row
    virtual void getrow(uint row, Array<uint>& columns, Array<real>& values) const;

    /// Set values for given row
    virtual void setrow(uint row, const Array<uint>& columns, const Array<real>& values);

    /// Set given rows to zero
    virtual void zero(uint m, const uint* rows);

    /// Set given rows to identity matrix
    virtual void ident(uint m, const uint* rows);

    // Matrix-vector product, y = Ax
    virtual void mult(const GenericVector& x, GenericVector& y, bool transposed=false) const;

    /// Multiply matrix by given number
    virtual const MTL4Matrix& operator*= (real a);

    /// Divide matrix by given number
    virtual const MTL4Matrix& operator/= (real a);

    /// Assignment operator
    virtual const GenericMatrix& operator= (const GenericMatrix& x)
    { error("Not implemented."); return *this; }

    //--- Special functions ---
    virtual LinearAlgebraFactory& factory() const;

    //--- Special MTL4 functions ---

    /// Return MTL4_sparse_matrix reference
    MTL4_sparse_matrix& mat() const;

    /// Assignment operator
    const MTL4Matrix& operator= (const MTL4Matrix& x)
    { error("Not implemented."); return *this; }

  private:

    // MTL4_FECrsMatrix pointer
    MTL4_sparse_matrix* A;
    
    // True if we don't own the matrix A points to
    bool is_view;

    mtl::matrix::inserter<MTL4_sparse_matrix, mtl::update_plus<real> >* ins;
  };

  LogStream& operator<< (LogStream& stream, const MTL4_sparse_matrix& A);

}

#endif
#endif
