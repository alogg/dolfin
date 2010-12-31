// Copyright (C) 2003-2009 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells, 2005-2010.
// Modified by Martin Sandve Alnes, 2008.
// Modified by Andre Massing, 2009.
//
// First added:  2003-11-28
// Last changed: 2010-12-28

#include <algorithm>
#include <map>
#include <utility>
#include <boost/assign/list_of.hpp>

#include <dolfin/adaptivity/Extrapolation.h>
#include <dolfin/common/utils.h>
#include <dolfin/common/Array.h>
#include <dolfin/common/NoDeleter.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/GenericDofMap.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/fem/UFC.h>
#include <dolfin/io/File.h>
#include <dolfin/la/GenericVector.h>
#include <dolfin/la/DefaultFactory.h>
#include <dolfin/log/log.h>
#include <dolfin/mesh/Vertex.h>
#include <dolfin/parameter/GlobalParameters.h>
#include "Expression.h"
#include "FunctionSpace.h"
#include "Function.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
Function::Function(const FunctionSpace& V)
  : _function_space(reference_to_no_delete_pointer(V)),
    allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Initialize vector
  init_vector();
}
//-----------------------------------------------------------------------------
Function::Function(boost::shared_ptr<const FunctionSpace> V)
  : _function_space(V),
    allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Initialize vector
  init_vector();
}
//-----------------------------------------------------------------------------
Function::Function(const FunctionSpace& V, GenericVector& x)
  : _function_space(reference_to_no_delete_pointer(V)),
    _vector(reference_to_no_delete_pointer(x)),
    allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Assertion uses '<=' to deal with sub-functions
  assert(V.dofmap().global_dimension() <= x.size());
}
//-----------------------------------------------------------------------------
Function::Function(boost::shared_ptr<const FunctionSpace> V,
                   boost::shared_ptr<GenericVector> x)
  : _function_space(V),
    _vector(x),
    allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Assertion uses '<=' to deal with sub-functions
  assert(V->dofmap().global_dimension() <= x->size());
}
//-----------------------------------------------------------------------------
Function::Function(boost::shared_ptr<const FunctionSpace> V,
                   GenericVector& x)
  : _function_space(V),
    _vector(reference_to_no_delete_pointer(x)),
    allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Assertion uses '<=' to deal with sub-functions
  assert(V->dofmap().global_dimension() <= x.size());
}
//-----------------------------------------------------------------------------
Function::Function(const FunctionSpace& V, std::string filename)
  : _function_space(reference_to_no_delete_pointer(V)),
    allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Initialize vector
  init_vector();

  // Read vector from file
  File file(filename);
  file >> *_vector;

  // Check size of vector
  if (_vector->size() != _function_space->dim())
    error("Unable to read Function from file, number of degrees of freedom (%d) does not match dimension of function space (%d).", _vector->size(), _function_space->dim());
}
//-----------------------------------------------------------------------------
Function::Function(boost::shared_ptr<const FunctionSpace> V,
                   std::string filename)
  : _function_space(V),
    allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Create vector
  DefaultFactory factory;
  _vector.reset(factory.create_vector());

  // Initialize vector
  init_vector();

  // Read vector from file
  File file(filename);
  file >> *_vector;

  // Check size of vector
  if (_vector->size() != _function_space->dim())
    error("Unable to read Function from file, number of degrees of freedom (%d) does not match dimension of function space (%d).", _vector->size(), _function_space->dim());
}
//-----------------------------------------------------------------------------
Function::Function(const Function& v)
  : allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Assign data
  *this = v;
}
//-----------------------------------------------------------------------------
Function::Function(const Function& v, uint i)
  : allow_extrapolation(dolfin::parameters["allow_extrapolation"])
{
  // Copy function space pointer
  this->_function_space = v[i]._function_space;

  // Copy vector pointer
  this->_vector = v[i]._vector;
}
//-----------------------------------------------------------------------------
Function::~Function()
{
  // Do nothing
}
//-----------------------------------------------------------------------------
const Function& Function::operator= (const Function& v)
{
  assert(v._vector);

  // Make a copy of all the data, or if v is a sub-function, then we collapse
  // the dof map and copy only the relevant entries from the vector of v.
  if (v._vector->size() == v._function_space->dim())
  {
    // Copy function space
    _function_space = v._function_space;

    // Copy vector
    _vector.reset(v._vector->copy());
  }
  else
  {
    // Create collapsed dof map
    std::map<uint, uint> collapsed_map;
    boost::shared_ptr<GenericDofMap> collapsed_dof_map(v._function_space->dofmap().collapse(collapsed_map, v._function_space->mesh()));

    // Create new FunctionsSpapce
    _function_space = v._function_space->collapse_sub_space(collapsed_dof_map);

    // FIXME: This assertion doesn't work in parallel
    //assert(collapsed_map.size() == _function_space->dofmap().global_dimension());
    //assert(collapsed_map.size() == _function_space->dofmap().local_dimension());

    // Create new vector (global)
    _vector.reset(v.vector().factory().create_vector());
    _vector->resize(collapsed_dof_map->global_dimension());

    // Get row indices of original and new vectors
    std::map<uint, uint>::const_iterator entry;
    std::vector<uint> new_rows(collapsed_map.size());
    Array<uint> old_rows(collapsed_map.size());
    uint i = 0;
    for (entry = collapsed_map.begin(); entry != collapsed_map.end(); ++entry)
    {
      new_rows[i]   = entry->first;
      old_rows[i++] = entry->second;
    }

    // Gather values into an Array
    Array<double> gathered_values;
    v.vector().gather(gathered_values, old_rows);

    // Set values in vector
    this->_vector->set(&gathered_values[0], collapsed_map.size(), &new_rows[0]);
    this->_vector->apply("insert");
  }

  return *this;
}
//-----------------------------------------------------------------------------
const Function& Function::operator= (const Expression& v)
{
  interpolate(v);
  return *this;
}
//-----------------------------------------------------------------------------
Function& Function::operator[] (uint i) const
{
  // Check if sub-Function is in the cache, otherwise create and add to cache
  boost::ptr_map<uint, Function>::iterator sub_function = sub_functions.find(i);
  if (sub_function != sub_functions.end())
    return *(sub_function->second);
  else
  {
    // Extract function subspace
    std::vector<uint> component = boost::assign::list_of(i);
    boost::shared_ptr<const FunctionSpace> sub_space(_function_space->extract_sub_space(component));

    // Insert sub-Function into map and return reference
    sub_functions.insert(i, new Function(sub_space, _vector));
    return *(sub_functions.find(i)->second);
  }
}
//-----------------------------------------------------------------------------
const FunctionSpace& Function::function_space() const
{
  assert(_function_space);
  return *_function_space;
}
//-----------------------------------------------------------------------------
boost::shared_ptr<const FunctionSpace> Function::function_space_ptr() const
{
  assert(_function_space);
  return _function_space;
}
//-----------------------------------------------------------------------------
GenericVector& Function::vector()
{
  // Check that this is not a sub function.
  if (_vector->size() != _function_space->dofmap().global_dimension())
  {
    cout << "Size of vector: " << _vector->size() << endl;
    cout << "Size of function space: " << _function_space->dofmap().global_dimension() << endl;
    error("You are attempting to access a non-const vector from a sub-Function.");
  }
  return *_vector;
}
//-----------------------------------------------------------------------------
const GenericVector& Function::vector() const
{
  assert(_vector);
  return *_vector;
}
//-----------------------------------------------------------------------------
bool Function::in(const FunctionSpace& V) const
{
  assert(_function_space);
  return _function_space.get() == &V;
}
//-----------------------------------------------------------------------------
dolfin::uint Function::geometric_dimension() const
{
  assert(_function_space);
  return _function_space->mesh().geometry().dim();
}
//-----------------------------------------------------------------------------
void Function::eval(Array<double>& values, const Array<double>& x) const
{
  assert(_function_space);

  // Find the cell that contains x
  const double* _x = x.data().get();
  const Point point(_function_space->mesh().geometry().dim(), _x);
  int id = _function_space->mesh().any_intersected_entity(point);

  // If not found, use the closest cell
  if (id == -1)
  {
    if (allow_extrapolation)
    {
      id = _function_space->mesh().closest_cell(point);
      cout << "Extrapolating function value at x = " << point << " (not inside domain)." << endl;
    }
    else
    {
      cout << "Evaluating at x = " << point << endl;
      error("Unable to evaluate function at given point (not inside domain). Set option 'allow_extrapolation' to allow extrapolation'.");
    }
  }

  // Create cell that contains point
  const Cell cell(_function_space->mesh(), id);
  const UFCCell ufc_cell(cell);

  // Call evaluate function
  eval(values, x, cell, ufc_cell);
}
//-----------------------------------------------------------------------------
void Function::eval(Array<double>& values,
                    const Array<double>& x,
                    const Cell& dolfin_cell,
                    const ufc::cell& ufc_cell) const
{
  // Developer note: work arrays/vectors are re-created each time this function
  //                 is called for thread-safety

  const FiniteElement& element = _function_space->element();

  // FIXME: Rather than computing num_tensor_entries, we could just use
  //        values.size()

  // Compute in tensor (one for scalar function, . . .)
  uint num_tensor_entries = 1;
  for (uint i = 0; i < element.value_rank(); i++)
    num_tensor_entries *= element.value_dimension(i);

  assert(values.size() == num_tensor_entries);

  // Create work vector for expansion coefficients
  std::vector<double> coefficients(element.space_dimension());

  // Restrict function to cell
  restrict(&coefficients[0], element, dolfin_cell, ufc_cell);

  // Create work vector for basis
  std::vector<double> basis(num_tensor_entries);

  // Initialise values
  for (uint j = 0; j < num_tensor_entries; ++j)
    values[j] = 0.0;

  // Compute linear combination
  for (uint i = 0; i < element.space_dimension(); ++i)
  {
    element.evaluate_basis(i, &basis[0], &x[0], ufc_cell);
    for (uint j = 0; j < num_tensor_entries; ++j)
      values[j] += coefficients[i]*basis[j];
  }
}
//-----------------------------------------------------------------------------
void Function::interpolate(const GenericFunction& v)
{
  // Gather off-process dofs
  v.gather();

  // Initialise vector
  init_vector();

  // Interpolate
  function_space().interpolate(*_vector, v);
}
//-----------------------------------------------------------------------------
void Function::extrapolate(const Function& v)
{
  Extrapolation::extrapolate(*this, v);
}
//-----------------------------------------------------------------------------
dolfin::uint Function::value_rank() const
{
  return _function_space->element().value_rank();
}
//-----------------------------------------------------------------------------
dolfin::uint Function::value_dimension(uint i) const
{
  return _function_space->element().value_dimension(i);
}
//-----------------------------------------------------------------------------
void Function::eval(Array<double>& values, const Array<double>& x,
                    const ufc::cell& ufc_cell) const
{
  assert(_function_space);

  // Check if UFC cell comes from mesh, otherwise redirect to point-based evaluation
  if (ufc_cell.mesh_identifier == (int) _function_space->mesh().id())
  {
    const Cell cell(_function_space->mesh(), ufc_cell.index);
    eval(values, x, cell, ufc_cell);
  }
  else
    eval(values, x);
}
//-----------------------------------------------------------------------------
void Function::restrict(double* w,
                        const FiniteElement& element,
                        const Cell& dolfin_cell,
                        const ufc::cell& ufc_cell) const
{
  assert(w);
  assert(_function_space);

  // Check if we are restricting to an element of this function space
  if (_function_space->has_element(element) && _function_space->has_cell(dolfin_cell))
  {
    // Get dofmap for cell
    const GenericDofMap& dofmap = _function_space->dofmap();
    const std::vector<uint>& dofs = dofmap.cell_dofs(dolfin_cell.index());

    // Pick values from vector(s)
    _vector->get_local(w, dofs.size(), &dofs[0]);
  }
  else
  {
    // Restrict as UFC function (by calling eval)
    restrict_as_ufc_function(w, element, dolfin_cell, ufc_cell);
  }
}
//-----------------------------------------------------------------------------
void Function::compute_vertex_values(Array<double>& vertex_values,
                                     const Mesh& mesh) const
{
  assert(&mesh == &_function_space->mesh());

  // Gather ghosts dofs
  gather();

  // Get finite element
  const FiniteElement& element = _function_space->element();

  // Local data for interpolation on each cell
  const uint num_cell_vertices = mesh.type().num_vertices(mesh.topology().dim());

  // Compute in tensor (one for scalar function, . . .)
  uint num_tensor_entries = 1;
  for (uint i = 0; i < element.value_rank(); i++)
    num_tensor_entries *= element.value_dimension(i);

  // Create vector to hold cell vertex values
  std::vector<double> cell_vertex_values(num_tensor_entries*num_cell_vertices);

  // Create vector for expansion coefficients
  std::vector<double> coefficients(element.space_dimension());

  // Interpolate vertex values on each cell (using last computed value if not
  // continuous, e.g. discontinuous Galerkin methods)
  UFCCell ufc_cell(mesh);
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    // Update to current cell
    ufc_cell.update(*cell);

    // Pick values from global vector
    restrict(&coefficients[0], element, *cell, ufc_cell);

    // Interpolate values at the vertices
    element.interpolate_vertex_values(&cell_vertex_values[0],
                                      &coefficients[0], ufc_cell);

    // Copy values to array of vertex values
    for (VertexIterator vertex(*cell); !vertex.end(); ++vertex)
    {
      for (uint i = 0; i < num_tensor_entries; ++i)
      {
        const uint local_index  = vertex.pos()*num_tensor_entries + i;
        const uint global_index = i*mesh.num_vertices() + vertex->index();
        vertex_values[global_index] = cell_vertex_values[local_index];
      }
    }
  }
}
//-----------------------------------------------------------------------------
void Function::gather() const
{
  if (MPI::num_processes() > 1)
    _vector->update_ghost_values();
}
//-----------------------------------------------------------------------------
void Function::init_vector()
{
  // Get global size
  const uint N = _function_space->dofmap().global_dimension();

  // Get range
  const std::pair<uint, uint> range = MPI::local_range(N);

  // Get local range
  const uint n0 = range.first;
  const uint n1 = range.second;
  const uint local_size = n1 - n0;

  // Get ghost vertices if vector is distributed
  std::vector<uint> ghost_indices;
  if (N != local_size)
    compute_ghost_indices(n0, n1, ghost_indices);

  // Create vector of dofs
  if (!_vector)
  {
    DefaultFactory factory;
    _vector.reset(factory.create_vector());
  }
  assert(_vector);

  // Initialize vector of dofs
  _vector->resize(N, local_size, ghost_indices);
  _vector->zero();
}
//-----------------------------------------------------------------------------
void Function::compute_ghost_indices(uint n0, uint n1,
                                     std::vector<uint>& ghost_indices) const
{
  // Clear data
  ghost_indices.clear();

  // Get mesh
  assert(_function_space);
  const Mesh& mesh = _function_space->mesh();

  // Get dof map
  const GenericDofMap& dofmap = _function_space->dofmap();

  // Dofs per cell and total dofs
  const uint num_dofs_per_cell = _function_space->element().space_dimension();

  // Iterate over local mesh and check which dofs are needed
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    // Get dofs on cell
    const std::vector<uint>& dofs = dofmap.cell_dofs(cell->index());

    for (uint d = 0; d < num_dofs_per_cell; ++d)
    {
      const uint dof = dofs[d];
      if (dof < n0 || dof >= n1)
      {
        // FIXME: Could we use dolfin::Set here? Or unordered_set?
        if (std::find(ghost_indices.begin(), ghost_indices.end(), dof) == ghost_indices.end())
          ghost_indices.push_back(dof);
      }
    }
  }
}
//-----------------------------------------------------------------------------
