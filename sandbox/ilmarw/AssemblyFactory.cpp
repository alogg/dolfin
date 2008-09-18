// Copyright (C) 2008 Ilmar Wilbers.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-05-19
// Last changed: 2008-05-19

#include "SparsityPattern.h"
#include "AssemblyMatrix.h"
//#include "AssemblyVector.h"
#include "AssemblyFactory.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
SparsityPattern* AssemblyFactory::create_pattern() const 
{
  return new SparsityPattern(); 
}
//-----------------------------------------------------------------------------
AssemblyMatrix* AssemblyFactory::create_matrix() const 
{ 
  AssemblyMatrix* pm = new AssemblyMatrix();
  return pm;
}
//-----------------------------------------------------------------------------
AssemblyVector* AssemblyFactory:: create_vector() const 
{ 
  return new AssemblyVector(); 
}
//-----------------------------------------------------------------------------

// Singleton instance
AssemblyFactory AssemblyFactory::factory;
