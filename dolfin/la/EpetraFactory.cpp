// Copyright (C) 2008 Martin Sandve Alnes, Kent-Andre Mardal and Johannes Ring.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2008-04-21

#ifdef HAS_TRILINOS

#include "EpetraSparsityPattern.h"
#include "EpetraMatrix.h"
#include "EpetraVector.h"
#include "EpetraFactory.h"

#include <Epetra_SerialComm.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
EpetraFactory:: EpetraFactory() {
  comm = new Epetra_SerialComm(); 
}
//-----------------------------------------------------------------------------
EpetraMatrix* EpetraFactory::createMatrix() const 
{ 
  return new EpetraMatrix();
}
//-----------------------------------------------------------------------------
EpetraSparsityPattern* EpetraFactory::createPattern() const 
{
  return new EpetraSparsityPattern(); 
}
//-----------------------------------------------------------------------------
EpetraVector* EpetraFactory::createVector() const 
{ 
  return new EpetraVector(); 
}
//-----------------------------------------------------------------------------
EpetraFactory EpetraFactory::epetrafactory;
//-----------------------------------------------------------------------------
Epetra_SerialComm& EpetraFactory::getSerialComm() { return *comm;}; 
//-----------------------------------------------------------------------------

#endif
