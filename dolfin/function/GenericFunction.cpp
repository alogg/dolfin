// Copyright (C) 2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2009-09-28
// Last changed: 2011-01-19

#include <string>
#include <dolfin/fem/FiniteElement.h>
#include "GenericFunction.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
GenericFunction::GenericFunction() : Variable("u", "a function")
{
  // Do nothing
}
//-----------------------------------------------------------------------------
GenericFunction::~GenericFunction()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
void GenericFunction::eval(Array<double>& values, const Array<double>& x,
                           const ufc::cell& cell) const
{
  // Redirect to simple eval
  eval(values, x);
}
//-----------------------------------------------------------------------------
void GenericFunction::eval(Array<double>& values, const Array<double>& x) const
{
  error("Missing eval() function (must be overloaded).");
}
//-----------------------------------------------------------------------------
double GenericFunction::operator() (double x)
{
  // Check that function is scalar
  if (value_rank() != 0)
    error("Function is not scalar.");

  // Set up Array arguments
  Array<double> value(1);
  Array<double> xx(1);
  xx[0] = x;

  // Call eval
  eval(value, xx);

  // Return value
  return value[0];
}
//-----------------------------------------------------------------------------
double GenericFunction::operator() (double x, double y)
{
  // Check that function is scalar
  if (value_rank() != 0)
    error("Function is not scalar.");

  // Set up Array arguments
  Array<double> value(1);
  Array<double> xx(2);
  xx[0] = x;
  xx[1] = y;

  // Call eval
  eval(value, xx);

  // Return value
  return value[0];
}
//-----------------------------------------------------------------------------
double GenericFunction::operator() (double x, double y, double z)
{
  // Check that function is scalar
  if (value_rank() != 0)
    error("Function is not scalar.");

  // Set up Array arguments
  Array<double> value(1);
  Array<double> xx(3);
  xx[0] = x;
  xx[1] = y;
  xx[2] = z;

  // Call eval
  eval(value, xx);

  // Return value
  return value[0];
}
//-----------------------------------------------------------------------------
double GenericFunction::operator() (const Point& p)
{
  return (*this)(p.x(), p.y(), p.z());
}
//-----------------------------------------------------------------------------
void GenericFunction::operator() (Array<double>& value,
                                  double x)
{
  // Set up Array argument
  Array<double> xx(1);
  xx[0] = x;

  // Call eval
  eval(value, xx);
}
//-----------------------------------------------------------------------------
void GenericFunction::operator() (Array<double>& value,
                                  double x, double y)
{
  // Set up Array argument
  Array<double> xx(2);
  xx[0] = x;
  xx[1] = y;

  // Call eval
  eval(value, xx);
}
//-----------------------------------------------------------------------------
void GenericFunction::operator() (Array<double>& value,
                                  double x, double y, double z)
{
  // Set up Array argument
  Array<double> xx(3);
  xx[0] = x;
  xx[1] = y;
  xx[2] = z;

  // Call eval
  eval(value, xx);
}
//-----------------------------------------------------------------------------
void GenericFunction::operator() (Array<double>& value, const Point& p)
{
  (*this)(value, p.x(), p.y(), p.z());
}
//-----------------------------------------------------------------------------
dolfin::uint GenericFunction::value_size() const
{
  uint size = 1;
  for (uint i = 0; i < value_rank(); ++i)
    size *= value_dimension(i);
  return size;
}
//-----------------------------------------------------------------------------
void GenericFunction::evaluate(double* values,
                               const double* coordinates,
                               const ufc::cell& cell) const
{
  assert(values);
  assert(coordinates);

  // Wrap data
  Array<double> _values(value_size(), values);
  const Array<double> x(cell.geometric_dimension, const_cast<double*>(coordinates));

  // Redirect to eval
  eval(_values, x, cell);
}
//-----------------------------------------------------------------------------
void GenericFunction::restrict_as_ufc_function(double* w,
                                               const FiniteElement& element,
                                               const Cell& dolfin_cell,
                                               const ufc::cell& ufc_cell) const
{
  assert(w);

  // Evaluate each dof to get the expansion coefficients
  for (uint i = 0; i < element.space_dimension(); ++i)
    w[i] = element.evaluate_dof(i, *this, ufc_cell);
}
//-----------------------------------------------------------------------------
