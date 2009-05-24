// Copyright (C) 2008 Martin Sandve Alnes, Kent-Andre Mardal and Johannes Ring.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2009.
//
// First added:  2008-04-21
// Last changed: 2009-05-23

#ifdef HAS_TRILINOS

#include <dolfin/common/types.h>
#include <dolfin/log/log.h>
#include "EpetraFactory.h"
#include "GenericSparsityPattern.h"
#include "EpetraSparsityPattern.h"

#include <Epetra_SerialComm.h>
#include <Epetra_FECrsGraph.h>

using namespace dolfin;
using dolfin::uint;

//-----------------------------------------------------------------------------
EpetraSparsityPattern::EpetraSparsityPattern() : rank(0), epetra_graph(0)
{
  dims[0] = 0;
  dims[1] = 0;
}
//-----------------------------------------------------------------------------
EpetraSparsityPattern::~EpetraSparsityPattern()
{
  delete epetra_graph;
}
//-----------------------------------------------------------------------------
void EpetraSparsityPattern::init(uint rank_, const uint* dims_)
{
  rank = rank_;

  if (rank == 1)
  {
    dims[0] = dims_[0];
  }
  else if (rank == 2)
  {
    dims[0] = dims_[0];
    dims[1] = dims_[1];

    EpetraFactory& f = EpetraFactory::instance();
    Epetra_SerialComm Comm = f.getSerialComm();

    Epetra_Map row_map(dims[0], 0, Comm);
    epetra_graph = new Epetra_FECrsGraph(Copy, row_map, 0);
  }
  else
    error("Illegal rank for Epetra sparsity pattern.");
}
//-----------------------------------------------------------------------------
void EpetraSparsityPattern::insert(const uint* num_rows,
                                   const uint * const * rows)
{
  if (rank == 2)
  {
    epetra_graph->InsertGlobalIndices(num_rows[0], reinterpret_cast<const int*>(rows[0]),
                                      num_rows[1], reinterpret_cast<const int*>(rows[1]));
  }
}
//-----------------------------------------------------------------------------
void EpetraSparsityPattern::sort()
{
  dolfin_not_implemented();
}
//-----------------------------------------------------------------------------
uint EpetraSparsityPattern::size(uint i) const
{
  if (rank == 1)
  {
    return dims[0];
  }
  
  if (rank == 2)
  {
    dolfin_assert(epetra_graph);
    if (i == 0)
      return epetra_graph->NumGlobalRows();
    else
      return epetra_graph->NumGlobalCols();
  }
  return 0;
}
//-----------------------------------------------------------------------------
uint EpetraSparsityPattern::num_nonzeros() const
{
  dolfin_not_implemented();
  return 0;
}
//-----------------------------------------------------------------------------
void EpetraSparsityPattern::num_nonzeros_diagonal(uint* num_nonzeros) const
{
  dolfin_not_implemented();
}
//-----------------------------------------------------------------------------
void EpetraSparsityPattern::num_nonzeros_off_diagonal(uint* num_nonzeros) const
{
  dolfin_not_implemented();
}
//-----------------------------------------------------------------------------
void EpetraSparsityPattern::apply()
{
  dolfin_assert(epetra_graph);

  // Could employ eg. OptimizeStorage. Not sure if this is wanted,
  // the graph would then depend on the equations, not only the method.

  EpetraFactory& f = EpetraFactory::instance();
  Epetra_SerialComm Comm = f.getSerialComm();

  Epetra_Map row_map(dims[0], 0, Comm);
  Epetra_Map col_map(dims[1], 0, Comm);

  epetra_graph->FillComplete (col_map, row_map);
}
//-----------------------------------------------------------------------------
Epetra_FECrsGraph& EpetraSparsityPattern:: pattern() const
{
  return *epetra_graph;
}
//-----------------------------------------------------------------------------

#endif
