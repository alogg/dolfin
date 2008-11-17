// Copyright (C) 2006-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Martin Sandve Alnes, 2008.
//
// First added:  2006-02-09
// Last changed: 2008-11-17

#include <dolfin/fem/FiniteElement.h>
#include <dolfin/log/dolfin_log.h>
#include <dolfin/mesh/Mesh.h>
#include "Constant.h"
#include "FunctionSpace.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Constant::Constant(double value)
  : Function(),
    _size(1), _values(0)
{
  _values = new double[1];
  _values[0] = value;
}
//-----------------------------------------------------------------------------
Constant::Constant(uint size, double value)
  : Function(),
    _size(size), _values(0)
{
  dolfin_assert(size > 0);

  _values = new double[size];
  for (uint i = 0; i < size; i++)
    _values[i] = value;
}
//-----------------------------------------------------------------------------
Constant::Constant(const std::vector<double>& values)
  : Function(),
    _size(values.size()), _values(0)
{
  dolfin_assert(values.size() > 0);

  _values = new double[values.size()];
  for (uint i = 0; i < values.size(); i++)
    _values[i] = values[i];
}
//-----------------------------------------------------------------------------
Constant::Constant(const std::vector<uint>& shape,
                   const std::vector<double>& values)
  : Function(),
    _size(0), _values(0)
{
  dolfin_assert(shape.size() > 0);
  dolfin_assert(values.size() > 0);

  // Compute size
  _size = 1;
  for (uint i = 0; i < shape.size(); i++)
    _size *= shape[i];

  // Copy values
  dolfin_assert(values.size() == _size);
  for (uint i = 0; i < values.size(); i++)
    _values[i] = values[i];
}
//-----------------------------------------------------------------------------
Constant::Constant(const Constant& c)
  : Function(),
    _size(0), _values(0)
{
  *this = c;
}
//-----------------------------------------------------------------------------
const Constant& Constant::operator= (const Constant& c)
{
  dolfin_assert(c._size > 0);
  dolfin_assert(c._values);

  delete _values;
  _size = c._size;
  _values = new double[c._size];
  for (uint i = 0; i < c._size; i++)
    _values[i] = c._values[i];

  return *this;
}
//-----------------------------------------------------------------------------
Constant::~Constant()
{
  delete [] _values;
}
//-----------------------------------------------------------------------------
void Constant::eval(double* values, const Data& data) const
{
  dolfin_assert(values);
  dolfin_assert(_values);

  // Copy values
  for (uint i = 0; i < _size; i++)
    values[i] = _values[i];
}
//-----------------------------------------------------------------------------
