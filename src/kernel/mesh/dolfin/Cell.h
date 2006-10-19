// Copyright (C) 2006 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2006-06-01
// Last changed: 2006-10-19

#ifndef __CELL_H
#define __CELL_H

#include <dolfin/MeshEntity.h>
#include <dolfin/MeshEntityIterator.h>

namespace dolfin
{

  /// A Cell is a MeshEntity of topological codimension 0.

  class Cell : public MeshEntity
  {
  public:

    /// Constructor
    Cell(Mesh& mesh, uint index) : MeshEntity(mesh, mesh.topology().dim(), index) {}

    /// Destructor
    ~Cell() {}
    
    /// Return alignment of given entity with respect to the cell
    inline uint alignment(uint dim, uint e) const { return _mesh.type().alignment(*this, dim, e); }

  };

  /// A CellIterator is a MeshEntityIterator of topological codimension 0.
  
  class CellIterator : public MeshEntityIterator
  {
  public:
    
    CellIterator(Mesh& mesh) : MeshEntityIterator(mesh, mesh.topology().dim()) {}
    CellIterator(MeshEntity& entity) : MeshEntityIterator(entity, entity.mesh().topology().dim()) {}
    CellIterator(MeshEntityIterator& it) : MeshEntityIterator(it, it->mesh().topology().dim()) {}

    inline Cell& operator*()
    { return static_cast<Cell&>(*static_cast<MeshEntityIterator>(*this)); }

    inline Cell* operator->()
    { return &static_cast<Cell&>(*static_cast<MeshEntityIterator>(*this)); }

  };    

}

#endif
