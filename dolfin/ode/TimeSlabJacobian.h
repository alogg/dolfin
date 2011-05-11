// Copyright (C) 2005-2008 Anders Logg
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN.  If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2005-01-28
// Last changed: 2008-04-08

#ifndef __TIME_SLAB_JACOBIAN_H
#define __TIME_SLAB_JACOBIAN_H

#include <dolfin/log/dolfin_log.h>
#include <dolfin/la/uBLASDenseMatrix.h>
#include <dolfin/la/uBLASVector.h>
#include <dolfin/la/uBLASKrylovMatrix.h>
#include <dolfin/la/Matrix.h>

namespace dolfin
{

  class ODE;
  class Method;
  class TimeSlab;

  /// This is the base class for Jacobians defined on mono- or
  /// multi-adaptive time slabs.

  class TimeSlabJacobian : public uBLASKrylovMatrix
  {
  public:

    /// Constructor
    TimeSlabJacobian(TimeSlab& timeslab);

    /// Destructor
    ~TimeSlabJacobian();

    /// Return number of rows (dim = 0) or columns (dim = 1)
    virtual uint size(uint dim) const = 0;

    /// Compute product y = Ax
    virtual void mult(const uBLASVector& x, uBLASVector& y) const = 0;

    /// (Re-)initialize computation of Jacobian
    virtual void init();

    /// Update dense copy of Jacobian
    void update();

    /// Return dense copy of Jacobian
    const uBLASDenseMatrix& matrix() const;

  protected:

    // The ODE
    ODE& ode;

    // Method, mcG(q) or mdG(q)
    const Method& method;

    // Dense copy of the Jacobian
    uBLASDenseMatrix A;

    // Vectors used to compute dense copy of the Jacobian
    uBLASVector ej, Aj;

  };

}

#endif
