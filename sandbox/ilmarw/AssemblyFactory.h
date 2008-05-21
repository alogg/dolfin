// Copyright (C) 2008 Ilmar Wilbers.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2007-11-30
// Last changed: 2007-12-06

#ifndef __PETSC_FACTORY_H
#define __PETSC_FACTORY_H

#include "AssemblyMatrix.h"
//#include "AssemblyVector.h"
#include "SparsityPattern.h"
#include "LinearAlgebraFactory.h"

namespace dolfin
{

  class AssemblyFactory : public LinearAlgebraFactory
  {
  public:

    /// Destructor
    virtual ~AssemblyFactory() {}

    /// Create empty matrix
    AssemblyMatrix* createMatrix() const;

    /// Create empty vector
    AssemblyVector* createVector() const;

    /// Create empty sparsity pattern 
    SparsityPattern* createPattern() const;

    /// Return singleton instance
    static AssemblyFactory& instance() 
    { return factory; }

  private:

    /// Private Constructor
    AssemblyFactory() {}
    static AssemblyFactory factory;

  };

}

#endif
