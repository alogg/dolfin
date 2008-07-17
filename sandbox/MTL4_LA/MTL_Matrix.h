// Copyright (C) 2008 Dag Lindbo
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-07-06
// Last changed: 2008-07-06

#ifndef __MTL_MATRIX_H
#define __MTL_MATRIX_H

#ifdef HAS_MTL

#include <boost/numeric/mtl/mtl.hpp>
//#include "/NOBACKUP/dag/mtl4/boost/numeric/mtl/mtl.hpp"

typedef mtl::compressed2D<double> MTL4_sparse_matrix;

namespace dolfin
{
  class MTL_Matrix
  {
  public:

    /// Create M x N matrix
    inline 
    MTL_Matrix(uint M, uint N, uint q) : 
      A(M,N), ins(0), nnz_row(q)
    {
      A = 0;
    }

    /// Copy constuctor
    MTL_Matrix(const MTL_Matrix& A);

    /// Destructor
    inline
    ~MTL_Matrix()
    {
      if(ins) delete ins;
    }

    /// Return "rank"
    inline
    uint rank(void) const
    {
      return 2;
    }

    /// Set all entries to zero and keep any sparse structure
    inline 
    void zero()
    {
      A *= 0;
    }

    /// Finalize assembly of tensor
    inline 
    void apply()
    {
      if(ins) delete ins;
      ins = 0;
    }

    /// Add block of values
    void add(const real* block, 
	     uint m, const uint* rows, 
	     uint n, const uint* cols)
    {
      if(!ins)
     	ins = new mtl::matrix::
     	  inserter<MTL4_sparse_matrix, mtl::update_plus<real> >(A, nnz_row);

      real val;
      for (uint i = 0; i < m; i++)
     	for (uint j = 0; j < n; j++)
     	  {
     	    val = block[i*n +j];
     	    if( val != 0 ) 
     	      (*ins)[rows[i]][cols[j]] <<  val;
     	  }
    }

    /// Add block of values
    void add_II(const real* block, 
		uint m, const uint* rows, 
		uint n, const uint* cols)
    {
      if(!ins)
    	ins = new mtl::matrix::
    	  inserter<MTL4_sparse_matrix, mtl::update_plus<real> >(A, nnz_row);

      error("Funcion add_II needs to be debugged.");
      //mtl::dense_vector<const uint> row(m,rows);
      //mtl::dense_vector<const uint> col(n,cols);
      //mtl::dense2D<const double> blk(m,n,block);

      //(*ins) << mtl::element_matrix(blk, row, col);
    }

    // does not compile yet...
    template <typename Block, typename Rows, typename Cols>
    void add(const Block& block, const Rows& rows, const Cols& cols, int a)
    {
      if(!ins)
	ins = new mtl::matrix::
	  inserter<MTL4_sparse_matrix, mtl::update_plus<real> >(A, nnz_row);

      (*ins) << mtl::element_matrix(block, rows, cols);
    }

    /// Add block of values
    inline 
    void add(const real* block, const uint* num_rows, 
	     const uint * const * rows)
    { 
      add(block, num_rows[0], rows[0], num_rows[1], rows[1]);
    }

    /// Return MTL__sparse_matrix reference
    inline
    const MTL4_sparse_matrix& mat() const
    {
      return A;
    }
    
    /// Return MTL__sparse_matrix reference
    inline
    MTL4_sparse_matrix& mat()
    {
      return A;
    }

  private:

    MTL4_sparse_matrix A;
    mtl::matrix::inserter<MTL4_sparse_matrix, mtl::update_plus<real> >* ins;
    uint nnz_row;
  };
}

#endif
#endif
