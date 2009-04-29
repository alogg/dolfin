// Copyright (C) 2005-2006 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2005-01-27
// Last changed: 2006-07-07

#ifndef __MULTI_ADAPTIVE_PRECONDITIONER_H
#define __MULTI_ADAPTIVE_PRECONDITIONER_H

#include <dolfin/la/uBLASPreconditioner.h>

namespace dolfin
{
  class ODE;
  class Method;
  class MultiAdaptiveTimeSlab;

  /// This class implements a preconditioner for the Newton system to
  /// be solved on a multi-adaptive time slab. The preconditioner just
  /// does simple forward propagation of values on internal elements
  /// of the time slab (without so much as looking at the Jacobian).

  class MultiAdaptivePreconditioner : public uBLASPreconditioner
  {
  public:

    /// Constructor
    MultiAdaptivePreconditioner(MultiAdaptiveTimeSlab& timeslab, const Method& method);

    /// Destructor
    ~MultiAdaptivePreconditioner();

    /// Solve linear system approximately for given right-hand side b
    void solve(uBLASVector& x, const uBLASVector& b) const;

  private:

    // The time slab
    MultiAdaptiveTimeSlab& ts;

    // Method, mcG(q) or mdG(q)
    const Method& method;

  };

}

#endif
