// Copyright (C) 2003 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#ifndef __SIMPLE_TIME_SLAB_H
#define __SIMPLE_TIME_SLAB_H

#include <dolfin/TimeSlab.h>

namespace dolfin {

  class Element;
  class TimeSlabData;
  class RHS;

  /// The simple version of the time slab.

  class SimpleTimeSlab : public TimeSlab {
  public:
    
    /// Create time slab, including one iteration
    SimpleTimeSlab(real t0, real t1, RHS& f, TimeSlabData& data);
    
    /// Destructor
    ~SimpleTimeSlab();
    
    /// Update time slab (iteration)
    void update(RHS& f, TimeSlabData& data);
    
  private:
    
    // Create new time slab
    void create(RHS& f, TimeSlabData& data);

    //--- Time slab data ---

    // List of elements within this time slab
    std::vector<Element*> elements;
    
  };

}

#endif
