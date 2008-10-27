// Copyright (C) 2007-2008 Garth N. Wells.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Anders Logg, 2008.
//
// First added:  2007-12-10
// Last changed: 2008-10-24

#include <ufc.h>
#include <dolfin/log/log.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/Function.h>
#include "Form.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Form::Form() : _ufc_form(0)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
Form::~Form()
{
  delete _ufc_form;
}
//-----------------------------------------------------------------------------
dolfin::uint Form::rank() const
{ 
  dolfin_assert(_ufc_form); 
  return _ufc_form->rank(); 
}
//-----------------------------------------------------------------------------
const Mesh& Form::mesh() const
{
  // Extract all meshes
  std::vector<const Mesh*> meshes;
  for (uint i = 0; i < _function_spaces.size(); i++)
    meshes.push_back(&_function_spaces[i]->mesh());
  for (uint i = 0; i < _coefficients.size(); i++)
    meshes.push_back(&_coefficients[i]->function_space().mesh());

  // Check that we have at least one mesh
  if (meshes.size() == 0)
    error("Unable to extract mesh from form (no mesh found).");

  // Check that all meshes are the same
  for (uint i = 1; i < meshes.size(); i++)
    if (meshes[i] != meshes[i - 1])
      error("Unable to extract mesh from form (nonmatching meshes for function spaces).");

  // Return first mesh
  dolfin_assert(meshes[0]);
  return *meshes[0];
}
//-----------------------------------------------------------------------------
const FunctionSpace& Form::function_space(uint i) const
{
  dolfin_assert(i < _function_spaces.size());
  return *_function_spaces[i];
}
//-----------------------------------------------------------------------------
const std::vector<const FunctionSpace*> Form::function_spaces() const
{
  std::vector<const FunctionSpace*> V;
  for (uint i = 0; i < _function_spaces.size(); ++i)
    V.push_back(_function_spaces[i].get());

  return V;
}
//-----------------------------------------------------------------------------
const Function& Form::coefficient(uint i) const
{
  dolfin_assert(i < _coefficients.size());
  return *_coefficients[i];
}
//-----------------------------------------------------------------------------
const std::vector<const Function*> Form::coefficients() const
{
  std::vector<const Function*> V;
  for (uint i = 0; i < _coefficients.size(); ++i)
    V.push_back(_coefficients[i].get());

  return V;
}
//-----------------------------------------------------------------------------
const ufc::form& Form::ufc_form() const
{
  dolfin_assert(_ufc_form);
  return *_ufc_form;
}
//-----------------------------------------------------------------------------
void Form::check() const
{
  // FIXME: Further checks are needed here, in particular checking that
  // FIXME: the element of each function space matches the corresponding
  // FIXME: element in the form (check signatures).

  // Check that the number of function spaces matches the rank of the form
  if (_ufc_form->rank() != _function_spaces.size())
    error("Form expects %d FunctionSpaces, only %d provided.",
          _ufc_form->rank(), _function_spaces.size());

  // Check that the number of coefficient functions matches the number expecte by the form
  if (_ufc_form->num_coefficients() != _coefficients.size())
    error("Form expects %d coefficient functions, only %d provided.",
          _ufc_form->num_coefficients(), _coefficients.size());
}
//-----------------------------------------------------------------------------
