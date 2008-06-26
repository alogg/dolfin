// Copyright (C) 2006 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-10-24
// Last changed: 2006-09-02

#include <dolfin/log/dolfin_log.h>
#include "NonlinearProblem.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
NonlinearProblem::NonlinearProblem()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
NonlinearProblem::~NonlinearProblem()
{
  // Do nothing 
}
//-----------------------------------------------------------------------------
void NonlinearProblem::form(GenericMatrix& A, GenericVector& b, const GenericVector& x)
{
  error("Nonlinear problem update for F(u) and J has not been supplied by user.");
}
//-----------------------------------------------------------------------------
void NonlinearProblem::F(GenericVector& b, const GenericVector& x)
{
  error("Nonlinear problem update for F(u) has not been supplied by user.");
}
//-----------------------------------------------------------------------------
void NonlinearProblem::J(GenericMatrix& A, const GenericVector& x)
{
  error("Nonlinear problem update for Jacobian has not been supplied by user.");
}
//-----------------------------------------------------------------------------
