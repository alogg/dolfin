// Copyright (C) 2003 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#ifndef __CELL_REF_DATA_H
#define __CELL_REF_DATA_H

#include <dolfin/Cell.h>

namespace dolfin {

  /// Cell refinement data
  class CellRefData {
  public:

    /// Create cell refinement data
    CellRefData() {
      marker = Cell::marked_according_to_ref;
    }
    
    /// The mark of the cell
    Cell::Marker marker;

  };

}

#endif
