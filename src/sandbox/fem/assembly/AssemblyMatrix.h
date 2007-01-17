// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-01-17
// Last changed: 2007-01-17

#ifndef __ASSEMBLY_MATRIX_H
#define __ASSEMBLY_MATRIX_H

#include <vector>
#include <map>

#include <dolfin/constants.h>

namespace dolfin
{

  /// Simple implementation of a GenericTensor for experimenting
  /// with new assembly. Not sure this will be used later but it
  /// might be useful.

  class AssemblyMatrix : public GenericTensor
  {
  public:

    /// Constructor
    AssemblyMatrix() : GenericTensor()
    {
      

    }

    /// Destructor
    ~AssemblyMatrix()
    {


    }

  private:

    std::vector<std::map<uint, real> > A;

  };

}

#endif
