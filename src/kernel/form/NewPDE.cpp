// Copyright (C) 2004 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#include <dolfin/NewPDE.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
NewPDE::NewPDE(BilinearForm& a, LinearForm& L) : a(a), L(L)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
NewPDE::~NewPDE()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
