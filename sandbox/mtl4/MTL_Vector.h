// Copyright (C) 2008 Dag Lindbo
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-07-06
// Last changed: 2008-07-09

#ifndef __MTL_VECTOR_H
#define __MTL_VECTOR_H

#ifdef HAS_MTL

#include <boost/numeric/mtl/mtl.hpp>
//#include "/NOBACKUP/dag/mtl4/boost/numeric/mtl/mtl.hpp"

typedef mtl::dense_vector<double> MTL4_vector;

namespace dolfin
{
  class MTL_Vector
  {
  public:

    /// Create vector
    inline 
    MTL_Vector(uint N) : v(N, 0.0)
    {}

    /// Copy constuctor
    MTL_Vector(const MTL_Vector& A);

    /// Destructor
    inline
    ~MTL_Vector(){}

    /// Return "rank"
    inline
    uint rank(void) const
    {
      return 1;
    }

    /// Set all entries to zero and keep any sparse structure
    inline 
    void zero()
    {
      v *= 0;
    }

    /// Finalize assembly of tensor
    inline 
    void apply(){}

    inline
    void add(const double* block, uint m, const uint* rows)
    {
      for (uint i = 0; i < m; i++)
	v[ rows[i] ] += block[i];
    }

    /// Add block of values
    inline 
    void add(const double* block, const uint* num_rows, 
	     const uint * const * rows)
    { 
      add(block, num_rows[0], rows[0]);
    }

    /// Return MTL_vectorx reference
    inline
    const MTL4_vector& vec() const
    {
      return v;
    }
    
    /// Return MTL_vector reference
    inline
    MTL4_vector& vec()
    {
      return v;
    }

  private:

    MTL4_vector v;
  };
}

#endif
#endif
