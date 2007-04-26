// Copyright (C) 2007 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2007-04-24
// Last changed: 2007-04-26

#include <dolfin/dolfin_log.h>
#include <dolfin/SubSystem.h>

using namespace dolfin;

//-----------------------------------------------------------------------------
SubSystem::SubSystem()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(uint sub_system)
{
  this->sub_system.push_back(sub_system);
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(uint sub_system, uint sub_sub_system)
{
  this->sub_system.push_back(sub_system);
  this->sub_system.push_back(sub_sub_system);
}
//-----------------------------------------------------------------------------
SubSystem::SubSystem(const Array<uint>& sub_system) : sub_system(sub_system)
{
  // Do nothing
}
//-----------------------------------------------------------------------------
const ufc::finite_element* SubSystem::extractFiniteElement
(const ufc::finite_element& finite_element) const
{
  return extractFiniteElement(finite_element, sub_system);
}
//-----------------------------------------------------------------------------
const ufc::dof_map* SubSystem::extractDofMap
(const ufc::dof_map& dof_map) const
{
  return extractDofMap(dof_map, sub_system);
}
//-----------------------------------------------------------------------------
const ufc::finite_element* SubSystem::extractFiniteElement
(const ufc::finite_element& finite_element, const Array<uint>& sub_system)
{
  // Check if there are any sub systems
  if (finite_element.num_sub_elements() == 0)
  {
    dolfin_error("Unable to extract sub system (there are no sub systems).");
  }

  // Check that a sub system has been specified
  if (sub_system.size() == 0)
  {
    dolfin_error("Unable to extract sub system (no sub system specified).");
  }
  
  // Check the number of available sub systems
  if (sub_system[0] >= finite_element.num_sub_elements())
  {
    dolfin_error2("Unable to extract sub system %d (only %d sub systems defined).",
                  sub_system[0], finite_element.num_sub_elements());
  }
  
  // Create sub system
  ufc::finite_element* sub_element = finite_element.create_sub_element(sub_system[0]);
  
  // Return sub system if sub sub system should not be extracted
  if (sub_system.size() == 1)
    return sub_element;

  // Otherwise, recursively extract the sub sub system
  Array<uint> sub_sub_system;
  for (uint i = 1; i < sub_system.size(); i++)
    sub_sub_system.push_back(sub_system[i]);
  const ufc::finite_element* sub_sub_element = extractFiniteElement(*sub_element, sub_sub_system);
  delete sub_element;

  return sub_sub_element;
}
//-----------------------------------------------------------------------------
const ufc::dof_map* SubSystem::extractDofMap
(const ufc::dof_map& dof_map, const Array<uint>& sub_system)
{
  if (sub_system.size() == 0)
    return &dof_map;

  return 0;
}
//-----------------------------------------------------------------------------
