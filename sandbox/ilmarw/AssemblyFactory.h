// Copyright (C) 2007 Ilmar Wilbers.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2008.
//
// First added:  2008-05-21
// Last changed: 2008-05-21

#ifndef __ASSEMBLY_FACTORY_H
#define __ASSEMBLY_FACTORY_H

#include "AssemblyMatrix.h"
#include "uBlasVector.h"
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
    AssemblyMatrix* createMatrix() const
    { return new AssemblyMatrix(); }

    /// Create empty vector
    uBlasVector* createVector() const
    { return new uBlasVector(); }

    /// Create empty sparsity pattern 
    SparsityPattern* createPattern() const
    { return new SparsityPattern(); }

    /// Return singleton instance
    static AssemblyFactory& instance() 
    { return factory; }

  private:

    /// Private Constructor
    AssemblyFactory() {}

    // Singleton instance
    static AssemblyFactory factory;
  };
}

// Initialise static data
AssemblyFactory AssemblyFactory::factory;

#endif
