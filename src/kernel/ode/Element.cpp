// Copyright (C) 2003 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#include <dolfin/dolfin_log.h>
#include <dolfin/RHS.h>
#include <dolfin/Element.h>

using namespace dolfin;

// Initialize static data
Vector dolfin::Element::f;

//-----------------------------------------------------------------------------
Element::Element(real t0, real t1, unsigned int q, unsigned int index) :
  t0(t0), t1(t1), q(q), _index(index)
{
  // Allocate the list of nodal values
  values = new real[q+1];

  for (unsigned int i = 0; i < q+1; i++)
    values[i] = 0;

  // Increase size of the common vector to the maximum required size
  // among all elements. The vector is short (of length max(q+1)), but
  // a lot of extra storage would be required if each element stored
  // its own vector.
  
  if ( (q+1) > f.size() )
    f.init(q+1);
}
//-----------------------------------------------------------------------------
Element::~Element()
{
  if ( values )
    delete [] values;
  values = 0;
}
//-----------------------------------------------------------------------------
real Element::value(unsigned int node) const
{
  dolfin_assert(node >= 0);
  dolfin_assert(node <= q);
  
  return values[node];
}
//-----------------------------------------------------------------------------
real Element::endval() const
{
  return values[q];
}
//-----------------------------------------------------------------------------
bool Element::within(real t) const
{
  return (t0 < t) && (t <= t1);
}
//-----------------------------------------------------------------------------
unsigned int Element::index() const
{
  return _index;
}
//-----------------------------------------------------------------------------
real Element::starttime() const
{
  return t0;
}
//-----------------------------------------------------------------------------
real Element::endtime() const
{
  return t1;
}
//-----------------------------------------------------------------------------
real Element::timestep() const
{
  return t1 - t0;
}
//-----------------------------------------------------------------------------
real Element::computeResidual(RHS& f)
{
  return dx() - f(_index, q, t1);
}
//-----------------------------------------------------------------------------
