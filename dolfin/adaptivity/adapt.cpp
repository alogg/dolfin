// Copyright (C) 2010-2011 Anders Logg, Marie Rognes and Garth N. Wells
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
// First added:  2010-02-10
// Last changed: 2012-01-20

#include <map>
#include <boost/shared_ptr.hpp>

#include <dolfin/common/types.h>

#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MeshEntity.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/refinement/LocalMeshRefinement.h>
#include <dolfin/refinement/UniformMeshRefinement.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/SubSpace.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/Function.h>
#include <dolfin/function/SpecialFacetFunction.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/DirichletBC.h>
#include <dolfin/fem/LinearVariationalProblem.h>
#include <dolfin/fem/NonlinearVariationalProblem.h>
#include <dolfin/plot/plot.h>
#include "ErrorControl.h"
#include "adapt.h"

using namespace dolfin;

// Common function for setting parent/child
template <typename T>
void set_parent_child(const T& parent, boost::shared_ptr<T> child)
{
  // Use a const_cast so we can set the parent/child
  T& _parent = const_cast<T&>(parent);

  // Set parent/child
  _parent.set_child(child);
  child->set_parent(reference_to_no_delete_pointer(_parent));
}
//-----------------------------------------------------------------------------
const Mesh& dolfin::adapt(const Mesh& mesh)
{
  // Skip refinement if already refined
  if (mesh.has_child())
  {
    dolfin_debug("Mesh has already been refined, returning child mesh.");
    return mesh.child();
  }

  // Refine uniformly
  boost::shared_ptr<Mesh> adapted_mesh(new Mesh());
  UniformMeshRefinement::refine(*adapted_mesh, mesh);

  // Initialize the entities initialized in mesh in adapted_mesh
  for (uint d = 0; d <= mesh.topology().dim(); ++d)
    if (mesh.num_entities(d) != 0)
      adapted_mesh->init(d);

  // Set parent / child
  set_parent_child(mesh, adapted_mesh);

  return *adapted_mesh;
}
//-----------------------------------------------------------------------------
const dolfin::Mesh& dolfin::adapt(const Mesh& mesh,
                                  const MeshFunction<bool>& cell_markers)
{
  // Skip refinement if already refined
  if (mesh.has_child())
  {
    dolfin_debug("Mesh has already been refined, returning child mesh.");
    return mesh.child();
  }

  // Call refinement algorithm
  boost::shared_ptr<Mesh> adapted_mesh(new Mesh());
  LocalMeshRefinement::refine(*adapted_mesh, mesh, cell_markers);

  // Initialize the entities initialized in mesh in adapted_mesh
  for (uint d = 0; d <= mesh.topology().dim(); ++d)
    if (mesh.num_entities(d) != 0)
      adapted_mesh->init(d);

  // Set parent / child
  set_parent_child(mesh, adapted_mesh);

  return *adapted_mesh;
}
//-----------------------------------------------------------------------------
const dolfin::FunctionSpace& dolfin::adapt(const FunctionSpace& space)
{
  dolfin_assert(space.mesh());

  // Refine mesh
  adapt(*space.mesh());

  // Refine space
  adapt(space, space.mesh()->child_shared_ptr());

  return space.child();
}
//-----------------------------------------------------------------------------
const dolfin::FunctionSpace& dolfin::adapt(const FunctionSpace& space,
                                           const MeshFunction<bool>& cell_markers)
{
  dolfin_assert(space.mesh());

  // Refine mesh
  adapt(*space.mesh(), cell_markers);

  // Refine space
  adapt(space, space.mesh()->child_shared_ptr());

  return space.child();
}
//-----------------------------------------------------------------------------
const dolfin::FunctionSpace& dolfin::adapt(const FunctionSpace& space,
                                           boost::shared_ptr<const Mesh> adapted_mesh)
{
  // Skip refinement if already refined
  if (space.has_child())
  {
    dolfin_debug("Function space has already been refined, returning child space.");
    return space.child();
  }

  // Create DOLFIN finite element and dofmap
  dolfin_assert(space.dofmap());
  dolfin_assert(space.element());
  boost::shared_ptr<const FiniteElement>
    refined_element(space.element()->create());
  boost::shared_ptr<const GenericDofMap>
    refined_dofmap(space.dofmap()->copy(*adapted_mesh));

  // Create new function space
  boost::shared_ptr<FunctionSpace>
    refined_space(new FunctionSpace(adapted_mesh, refined_element, refined_dofmap));

  // Set parent / child
  set_parent_child(space, refined_space);

  return *refined_space;
}
//-----------------------------------------------------------------------------
const dolfin::Function& dolfin::adapt(const Function& function,
                                      boost::shared_ptr<const Mesh> adapted_mesh,
                                      bool interpolate)
{
  // Skip refinement if already refined
  if (function.has_child())
  {
    dolfin_debug("Function has already been refined, returning child function.");
    return function.child();
  }

  // Refine function space
  boost::shared_ptr<const FunctionSpace> space = function.function_space();
  adapt(*space, adapted_mesh);
  boost::shared_ptr<const FunctionSpace>
    refined_space = space->child_shared_ptr();

  // Create new function on refined space and interpolate
  boost::shared_ptr<Function> refined_function(new Function(refined_space));
  if (interpolate)
    refined_function->interpolate(function);

  // Set parent / child
  set_parent_child(function, refined_function);

  return *refined_function;
}
//-----------------------------------------------------------------------------
const dolfin::GenericFunction& dolfin::adapt(const GenericFunction& function,
                                             boost::shared_ptr<const Mesh> adapted_mesh)
{
  // Try casting to a function
  const Function* f = dynamic_cast<const Function*>(&function);
  if (f)
    return adapt(*f, adapted_mesh);
  else
    return function;
}
//-----------------------------------------------------------------------------
const dolfin::Form& dolfin::adapt(const Form& form,
                                  boost::shared_ptr<const Mesh> adapted_mesh,
                                  bool adapt_coefficients)
{
  // Skip refinement if already refined
  if (form.has_child())
  {
    dolfin_debug("Form has already been refined, returning child form.");
    return form.child();
  }

  // Get data
  std::vector<boost::shared_ptr<const FunctionSpace> >
    spaces = form.function_spaces();
  std::vector<boost::shared_ptr<const GenericFunction> >
    coefficients = form.coefficients();
  boost::shared_ptr<const ufc::form> ufc_form = form.ufc_form();

  // Refine function spaces
  std::vector<boost::shared_ptr<const FunctionSpace> > refined_spaces;
  for (uint i = 0; i < spaces.size(); i++)
  {
    const FunctionSpace& space = *spaces[i];
    adapt(space, adapted_mesh);
    refined_spaces.push_back(space.child_shared_ptr());
  }

  // Refine coefficients:
  std::vector<boost::shared_ptr<const GenericFunction> > refined_coefficients;
  for (uint i = 0; i < coefficients.size(); i++)
  {
    // Try casting to Function
    const Function*
      function = dynamic_cast<const Function*>(coefficients[i].get());

    if (function)
    {
      adapt(*function, adapted_mesh, adapt_coefficients);
      refined_coefficients.push_back(function->child_shared_ptr());
    } else
      refined_coefficients.push_back(coefficients[i]);
  }

  /// Create new form (constructor used from Python interface)
  boost::shared_ptr<Form> refined_form(new Form(ufc_form,
                                                refined_spaces,
                                                refined_coefficients));

  /// Attach mesh
  refined_form->set_mesh(adapted_mesh);

  // Attached refined sub domains
  const MeshFunction<uint>* cell_domains = form.cell_domains_shared_ptr().get();
  if (cell_domains)
  {
    adapt(*cell_domains, adapted_mesh);
    refined_form->dx = cell_domains->child_shared_ptr();
  }
  const MeshFunction<uint>* exterior_domains = form.exterior_facet_domains_shared_ptr().get();
  if (exterior_domains)
  {
    adapt(*exterior_domains, adapted_mesh);
    refined_form->ds = exterior_domains->child_shared_ptr();
  }
  const MeshFunction<uint>* interior_domains = form.interior_facet_domains_shared_ptr().get();
  if (interior_domains)
  {
    adapt(*interior_domains, adapted_mesh);
    refined_form->dS = interior_domains->child_shared_ptr();
  }

  // Set parent / child
  set_parent_child(form, refined_form);

  return *refined_form;
}
//-----------------------------------------------------------------------------
const dolfin::LinearVariationalProblem&
dolfin::adapt(const LinearVariationalProblem& problem,
              boost::shared_ptr<const Mesh> adapted_mesh)
{
  // Skip refinement if already refined
  if (problem.has_child())
  {
    dolfin_debug("Linear variational problem has already been refined, returning child problem.");
    return problem.child();
  }

  // Get data
  boost::shared_ptr<const Form> a = problem.bilinear_form();
  boost::shared_ptr<const Form> L = problem.linear_form();
  boost::shared_ptr<const Function> u = problem.solution();
  std::vector<boost::shared_ptr<const BoundaryCondition> > bcs = problem.bcs();

  // Refine forms
  dolfin_assert(a);
  dolfin_assert(L);
  adapt(*a, adapted_mesh);
  adapt(*L, adapted_mesh);

  // FIXME: Note const-cast here, don't know how to get around it

  // Refine solution variable
  dolfin_assert(u);
  adapt(*u, adapted_mesh);
  boost::shared_ptr<Function> refined_u =
    reference_to_no_delete_pointer(const_cast<Function&>(u->child()));

  // Refine bcs
  boost::shared_ptr<const FunctionSpace> V(problem.trial_space());
  std::vector<boost::shared_ptr<const BoundaryCondition> > refined_bcs;
  for (uint i = 0; i < bcs.size(); i++)
  {
    const DirichletBC* bc = dynamic_cast<const DirichletBC*>(bcs[i].get());
    if (bc != 0)
    {
      dolfin_assert(V);
      adapt(*bc, adapted_mesh, *V);
      refined_bcs.push_back(bc->child_shared_ptr());
    }
    else
    {
      dolfin_error("adapt.cpp",
                   "adapt linear variational problem",
                   "Only implemented for Dirichlet boundary conditions");
    }
  }

  // Create new problem
  dolfin_assert(a);
  dolfin_assert(L);
  dolfin_assert(u);
  boost::shared_ptr<LinearVariationalProblem>
    refined_problem(new LinearVariationalProblem(a->child_shared_ptr(),
                                                 L->child_shared_ptr(),
                                                 refined_u,
                                                 refined_bcs));

  // Set parent / child
  set_parent_child(problem, refined_problem);

  return *refined_problem;
}
//-----------------------------------------------------------------------------
const dolfin::NonlinearVariationalProblem&
dolfin::adapt(const NonlinearVariationalProblem& problem,
              boost::shared_ptr<const Mesh> adapted_mesh)
{
  // Skip refinement if already refined
  if (problem.has_child())
  {
    dolfin_debug("Nonlinear variational problem has already been refined, returning child problem.");
    return problem.child();
  }

  // Get data
  boost::shared_ptr<const Form> F = problem.residual_form();
  boost::shared_ptr<const Form> J = problem.jacobian_form();
  boost::shared_ptr<const Function> u = problem.solution();
  std::vector<boost::shared_ptr<const BoundaryCondition> > bcs = problem.bcs();

  // Refine forms
  dolfin_assert(F);
  adapt(*F, adapted_mesh);
  if (J)
    adapt(*J, adapted_mesh);

  // FIXME: Note const-cast here, don't know how to get around it

  // Refine solution variable
  dolfin_assert(u);
  adapt(*u, adapted_mesh);
  boost::shared_ptr<Function> refined_u =
    reference_to_no_delete_pointer(const_cast<Function&>(u->child()));

  // Refine bcs
  boost::shared_ptr<const FunctionSpace> V(problem.trial_space());
  std::vector<boost::shared_ptr<const BoundaryCondition> > refined_bcs;
  for (uint i = 0; i < bcs.size(); i++)
  {
    const DirichletBC* bc = dynamic_cast<const DirichletBC*>(bcs[i].get());
    if (bc != 0)
    {
      dolfin_assert(V);
      adapt(*bc, adapted_mesh, *V);
      refined_bcs.push_back(bc->child_shared_ptr());
    }
    else
    {
      dolfin_error("adapt.cpp",
                   "adapt nonlinear variational problem",
                   "Only implemented for Dirichlet boundary conditions");
    }
  }

  // Create new problem
  dolfin_assert(F);
  dolfin_assert(u);
  boost::shared_ptr<NonlinearVariationalProblem> refined_problem;
  if (J)
    refined_problem.reset(new NonlinearVariationalProblem(F->child_shared_ptr(),
                                                          refined_u,
                                                          refined_bcs,
                                                          J->child_shared_ptr()));
  else
    refined_problem.reset(new NonlinearVariationalProblem(F->child_shared_ptr(),
                                                          refined_u,
                                                          refined_bcs));

  // Set parent / child
  set_parent_child(problem, refined_problem);

  return *refined_problem;
}
//-----------------------------------------------------------------------------
const dolfin::DirichletBC& dolfin::adapt(const DirichletBC& bc,
                                    boost::shared_ptr<const Mesh> adapted_mesh,
                                    const FunctionSpace& S)
{
  dolfin_assert(adapted_mesh);

  // Skip refinement if already refined
  if (bc.has_child())
  {
    dolfin_debug("DirichletBC has already been refined, returning child.");
    return bc.child();
  }

  boost::shared_ptr<const FunctionSpace> W = bc.function_space();
  dolfin_assert(W);

  // Refine function space
  const std::vector<uint> component = W->component();
  boost::shared_ptr<const FunctionSpace> V;
  if (component.size() == 0)
  {
    adapt(*W, adapted_mesh);
    V = W->child_shared_ptr();
  }
  else
  {
    adapt(S, adapted_mesh);
    V.reset(new SubSpace(S.child(), component));
  }

  // Get refined value
  const GenericFunction& g = adapt(*(bc.value()), adapted_mesh);
  boost::shared_ptr<const GenericFunction> g_ptr(reference_to_no_delete_pointer(g));

  // Extract user_sub_domain
  boost::shared_ptr<const SubDomain> user_sub_domain = bc.user_sub_domain();

  // Create refined boundary condition
  boost::shared_ptr<DirichletBC> refined_bc;
  if (user_sub_domain)
  {
    // Use user defined sub domain if defined
    refined_bc.reset(new DirichletBC(V, g_ptr, user_sub_domain, bc.method()));
  }
  else
  {
    // Extract markers
    const std::vector<std::pair<uint, uint> >& markers = bc.markers();

    // Create refined markers
    dolfin_assert(W->mesh());
    std::vector<std::pair<uint, uint> > refined_markers;
    adapt_markers(refined_markers, *adapted_mesh, markers, *W->mesh());

    refined_bc.reset(new DirichletBC(V, g_ptr, refined_markers, bc.method()));
  }

  // Set parent / child
  set_parent_child(bc, refined_bc);

  return *refined_bc;
}
//-----------------------------------------------------------------------------
const dolfin::ErrorControl& dolfin::adapt(const ErrorControl& ec,
                                    boost::shared_ptr<const Mesh> adapted_mesh,
                                    bool adapt_coefficients)
{
  dolfin_assert(adapted_mesh);

  // Skip refinement if already refined
  if (ec.has_child())
  {
    dolfin_debug("ErrorControl has already been refined, returning child problem.");
    return ec.child();
  }

  // Refine data
  adapt(*ec._residual, adapted_mesh, adapt_coefficients);
  adapt(*ec._L_star, adapted_mesh, adapt_coefficients);
  adapt(*ec._a_star, adapted_mesh, adapt_coefficients);
  adapt(*ec._a_R_T, adapted_mesh, adapt_coefficients);
  adapt(*ec._L_R_T, adapted_mesh, adapt_coefficients);
  adapt(*ec._a_R_dT, adapted_mesh, adapt_coefficients);
  adapt(*ec._L_R_dT, adapted_mesh, adapt_coefficients);
  adapt(*ec._eta_T, adapted_mesh, adapt_coefficients);

  // Create refined error control
  boost::shared_ptr<ErrorControl>
    refined_ec(new ErrorControl(ec._a_star->child_shared_ptr(),
                                ec._L_star->child_shared_ptr(),
                                ec._residual->child_shared_ptr(),
                                ec._a_R_T->child_shared_ptr(),
                                ec._L_R_T->child_shared_ptr(),
                                ec._a_R_dT->child_shared_ptr(),
                                ec._L_R_dT->child_shared_ptr(),
                                ec._eta_T->child_shared_ptr(),
                                ec._is_linear));

  // Set parent / child
  set_parent_child(ec, refined_ec);

  return *refined_ec;
}
//-----------------------------------------------------------------------------
const dolfin::MeshFunction<dolfin::uint>&
dolfin::adapt(const MeshFunction<uint>& mesh_function,
              boost::shared_ptr<const Mesh> adapted_mesh)
{
  // Skip refinement if already refined
  if (mesh_function.has_child())
  {
    dolfin_debug("MeshFunction has already been refined, returning child");
    return mesh_function.child();
  }

  const Mesh& mesh = mesh_function.mesh();
  const uint dim = mesh.topology().dim();

  // Extract parent map from data of refined mesh
  boost::shared_ptr<MeshFunction<unsigned int> > parent;
  if (mesh_function.dim() == dim)
    parent = adapted_mesh->data().mesh_function("parent_cell");
  else if (mesh_function.dim() == (dim - 1))
    parent = adapted_mesh->data().mesh_function("parent_facet");
  else
    dolfin_not_implemented();

  // Check that parent map exists
  if (!parent.get())
  {
    dolfin_error("adapt.cpp",
                 "adapt mesh function",
                 "Unable to extract information about parent mesh entites");
  }

  // Use very large value as 'undefined'
  const uint undefined = std::numeric_limits<unsigned int>::max();

  // Map values of mesh function into refined mesh function
  boost::shared_ptr<MeshFunction<uint> >
    adapted_mesh_function(new MeshFunction<uint>(*adapted_mesh,
                                                 mesh_function.dim()));
  for (uint i = 0; i < adapted_mesh_function->size(); i++)
  {
    const uint parent_index = (*parent)[i];
    if (parent_index < mesh_function.size())
      (*adapted_mesh_function)[i] = mesh_function[parent_index];
    else
      (*adapted_mesh_function)[i] = undefined;
  }

  // Set parent / child relations
  set_parent_child(mesh_function, adapted_mesh_function);

  // Return refined mesh function
  return *adapted_mesh_function;
}
//-----------------------------------------------------------------------------
void dolfin::adapt_markers(std::vector<std::pair<uint, uint> >& refined_markers,
                           const Mesh& adapted_mesh,
                           const std::vector<std::pair<uint, uint> >& markers,
                           const Mesh& mesh)
{

  // Extract parent map from data of refined mesh
  boost::shared_ptr<MeshFunction<unsigned int> > parent_cells = \
    adapted_mesh.data().mesh_function("parent_cell");
  boost::shared_ptr<MeshFunction<unsigned int> > parent_facets = \
    adapted_mesh.data().mesh_function("parent_facet");

  // Check that parent maps exist
  if (!parent_cells.get() || !parent_facets.get())
  {
    dolfin_error("adapt.cpp",
                 "adapt markers",
                 "Unable to extract information about parent mesh entites");
  }

  // Create map (parent_cell, parent_local_facet) -> [(child_cell,
  // child_local_facet), ...] for boundary facets
  std::pair<uint, uint> child;
  std::pair<uint, uint> parent;
  std::map< std::pair<uint, uint>,
    std::vector< std::pair<uint, uint> > > children;

  const uint D = mesh.topology().dim();
  for (FacetIterator facet(adapted_mesh); !facet.end(); ++facet)
  {
    // Ignore interior facets
    if (facet->num_entities(D) == 2)
      continue;

    // Extract cell and local facet number
    Cell cell(adapted_mesh, facet->entities(D)[0]);
    const uint local_facet = cell.index(*facet);

    child.first = cell.index();
    child.second = local_facet;

    // Extract parent cell
    Cell parent_cell(mesh, (*parent_cells)[cell]);

    // Extract (global) index of parent facet
    const uint parent_facet_index = (*parent_facets)[*facet];

    // Extract local number of parent facet wrt parent cell
    Facet parent_facet(mesh, parent_facet_index);
    const uint parent_local_facet = parent_cell.index(parent_facet);

    parent.first = parent_cell.index();
    parent.second = parent_local_facet;

    // Add this (cell, local_facet) to list of child facets
    children[parent].push_back(child);

  }

  // Use above map to construct refined markers
  std::vector<std::pair<uint, uint> >  child_facets;
  std::vector<std::pair<uint, uint> >::const_iterator it;
  for (it = markers.begin(); it != markers.end(); ++it)
  {
    child_facets = children[*it];
    for (uint k = 0; k < child_facets.size(); k++)
    {
      refined_markers.push_back(child_facets[k]);
    }
  }
}
//-----------------------------------------------------------------------------
