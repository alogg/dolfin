// Copyright (C) 2012 Anders Logg
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2012-04-11
// Last changed: 2012-04-13

#ifndef __CSG_PRIMITIVE_H
#define __CSG_PRIMITIVE_H

#include "CSGGeometry.h"

namespace dolfin
{

  /// Base class for Constructive Solid Geometry (CSG) primitives.

  class CSGPrimitive : public CSGGeometry
  {
  public:
    enum PrimitiveType { Box, Sphere, Cone, Circle, Rectangle, Polygon };
    virtual PrimitiveType getPrimitiveType() = 0;
  };

}

#endif
