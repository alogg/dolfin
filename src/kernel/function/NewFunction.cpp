// Copyright (C) 2005 Johan Hoffman and Anders Logg.
// Licensed under the GNU GPL Version 2.

#include <dolfin/Mesh.h>
#include <dolfin/NewFiniteElement.h>
#include <dolfin/NewVector.h>
#include <dolfin/NewFunction.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
NewFunction::NewFunction() : data(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
NewFunction::NewFunction(const Mesh& mesh, const NewFiniteElement& element,
			 NewVector& x) : data(0)
{
  // Create function data
  data = new Data(mesh, element, x);
  
  // Set name and label
  rename("u", "An unspecified function");
}
//-----------------------------------------------------------------------------
NewFunction::~NewFunction()
{
  if ( data ) delete data;
}
//-----------------------------------------------------------------------------
void NewFunction::project(const Cell& cell, const NewFiniteElement& element,
			  real c[]) const
{
  // Check if function is user-defined
  if ( !data )
  {
    for (uint i = 0; i < element.spacedim(); i++)
      c[i] = (*this)(element.coord(i, cell, data->mesh));

    return;
  }

  // Check if we're computing the projection onto a cell of the same
  // mesh for the same element (the easy case...)
  if ( &(cell.mesh()) == &(data->mesh) && &element == &(data->element) )
  {
    // FIXME: Assumes uniprocessor case. Why isn't there a function
    // FIXME: VecGetValues() in PETSc? Possible answer: since if we're
    // FIXME: doing this in parallel we only want to access this
    // FIXME: processor's data anyway.

    // FIXME: If we know that the values are stored element-by-element
    // FIXME: in x, then we can optimize by just calling
    // FIXME: element::dof() one time with i = 0.

    real *values = data->x.array();
    for (uint i = 0; i < element.spacedim(); i++)
      c[i] = values[element.dof(i, cell, data->mesh)];
    data->x.restore(values);

    return;
  }

  // Need to compute projection between different spaces
  dolfin_error("Projection between different finite element spaces not implemented.");
}
//-----------------------------------------------------------------------------
real NewFunction::operator()(const Point& p) const
{
  dolfin_error("User-defined function evaluation not implemented.");
  return 0.0;
}
//-----------------------------------------------------------------------------
