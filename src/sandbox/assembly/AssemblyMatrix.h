// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-17
// Last changed: 2007-03-01

#ifndef __ASSEMBLY_MATRIX_H
#define __ASSEMBLY_MATRIX_H

#include <vector>
#include <map>

#include <dolfin/constants.h>
#include <dolfin/dolfin_log.h>

namespace dolfin
{

  /// Simple implementation of a GenericTensor for experimenting
  /// with new assembly. Not sure this will be used later but it
  /// might be useful.

  class AssemblyMatrix : public GenericTensor
  {
  public:

    /// Constructor
    AssemblyMatrix() : GenericTensor(), dims(0)
    {
      dims = new uint[2];
    }

    /// Destructor
    ~AssemblyMatrix()
    {
      delete [] dims;
    }

    /// Initialize zero tensor of given rank and dimensions
    void init(uint rank, uint* dims)
    {
      // Check that the rank is 2
      if ( rank != 2 )
        dolfin_error1("Illegal tensor rank (%d) for matrix. Rank must be 2.", rank);

      // Initialize matrix
      init(dims[0], dims[1]);

      // Save dimensions
      this->dims[0] = dims[0];
      this->dims[1] = dims[1];
    }
    
    /// Return size of given dimension
    virtual uint size(const uint dim) const
    {
      return dims[dim];
    }

    /// Initialize M x N matrix
    void init(uint M, uint N)
    {
      // Set number of rows
      A.resize(M);
      
      // Initialize with zeros
      for (uint i = 0; i < M; i++)
        for (std::map<uint, real>::iterator it = A[i].begin(); it != A[i].end(); it++)
          it->second = 0.0;
    }

  private:

    // The matrix representation
    std::vector<std::map<uint, real> > A;

    // The size of the matrix
    uint* dims;

  };

}

#endif
