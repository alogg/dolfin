// Copyright (C) 2011 Anders Logg
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
// First added:  2011-06-22
// Last changed: 2011-07-17

#include <dolfin/common/NoDeleter.h>
#include <dolfin/function/Function.h>
#include "Form.h"
#include "NonlinearVariationalProblem.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
NonlinearVariationalProblem::
NonlinearVariationalProblem(const Form& F,
                            Function& u)
  : Hierarchical<NonlinearVariationalProblem>(*this),
    _F(reference_to_no_delete_pointer(F)),
    _u(reference_to_no_delete_pointer(u))
{
  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
NonlinearVariationalProblem::
NonlinearVariationalProblem(const Form& F,
                            Function& u,
                            const Form& J)
  : Hierarchical<NonlinearVariationalProblem>(*this),
    _F(reference_to_no_delete_pointer(F)),
    _J(reference_to_no_delete_pointer(J)),
    _u(reference_to_no_delete_pointer(u))
{
  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
NonlinearVariationalProblem::
NonlinearVariationalProblem(const Form& F,
                            Function& u,
                            const BoundaryCondition& bc)
  : Hierarchical<NonlinearVariationalProblem>(*this),
    _F(reference_to_no_delete_pointer(F)),
    _u(reference_to_no_delete_pointer(u))
{
  // Store boundary condition
  _bcs.push_back(reference_to_no_delete_pointer(bc));

  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
NonlinearVariationalProblem::
NonlinearVariationalProblem(const Form& F,
                            Function& u,
                            const BoundaryCondition& bc,
                            const Form& J)
  : Hierarchical<NonlinearVariationalProblem>(*this),
    _F(reference_to_no_delete_pointer(F)),
    _J(reference_to_no_delete_pointer(J)),
    _u(reference_to_no_delete_pointer(u))
{
  // Store boundary condition
  _bcs.push_back(reference_to_no_delete_pointer(bc));

  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
NonlinearVariationalProblem::
NonlinearVariationalProblem(const Form& F,
                            Function& u,
                            std::vector<const BoundaryCondition*> bcs)
  : Hierarchical<NonlinearVariationalProblem>(*this),
    _F(reference_to_no_delete_pointer(F)),
    _u(reference_to_no_delete_pointer(u))
{
  // Store boundary conditions
  for (std::size_t i = 0; i < bcs.size(); ++i)
    _bcs.push_back(reference_to_no_delete_pointer(*bcs[i]));

  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
NonlinearVariationalProblem::
NonlinearVariationalProblem(const Form& F,
                            Function& u,
                            std::vector<const BoundaryCondition*> bcs,
                            const Form& J)
  : Hierarchical<NonlinearVariationalProblem>(*this),
    _F(reference_to_no_delete_pointer(F)),
    _J(reference_to_no_delete_pointer(J)),
    _u(reference_to_no_delete_pointer(u))
{
  // Store boundary conditions
  for (std::size_t i = 0; i < bcs.size(); ++i)
    _bcs.push_back(reference_to_no_delete_pointer(*bcs[i]));

  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
NonlinearVariationalProblem::
NonlinearVariationalProblem(boost::shared_ptr<const Form> F,
                            boost::shared_ptr<Function> u,
                            std::vector<boost::shared_ptr<const BoundaryCondition> > bcs)
  : Hierarchical<NonlinearVariationalProblem>(*this),
    _F(F), _u(u)
{
  // Store boundary conditions
  for (std::size_t i = 0; i < bcs.size(); ++i)
    _bcs.push_back(bcs[i]);

  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
NonlinearVariationalProblem::
NonlinearVariationalProblem(boost::shared_ptr<const Form> F,
                            boost::shared_ptr<Function> u,
                            std::vector<boost::shared_ptr<const BoundaryCondition> > bcs,
                            boost::shared_ptr<const Form> J)
  : Hierarchical<NonlinearVariationalProblem>(*this),
    _F(F), _J(J), _u(u)
{
  // Store boundary conditions
  for (std::size_t i = 0; i < bcs.size(); ++i)
    _bcs.push_back(bcs[i]);

  // Check forms
  check_forms();
}
//-----------------------------------------------------------------------------
boost::shared_ptr<const Form> NonlinearVariationalProblem::residual_form() const
{
  return _F;
}
//-----------------------------------------------------------------------------
boost::shared_ptr<const Form> NonlinearVariationalProblem::jacobian_form() const
{
  return _J;
}
//-----------------------------------------------------------------------------
boost::shared_ptr<Function> NonlinearVariationalProblem::solution()
{
  return _u;
}
//-----------------------------------------------------------------------------
boost::shared_ptr<const Function> NonlinearVariationalProblem::solution() const
{
  return _u;
}
//-----------------------------------------------------------------------------
std::vector<boost::shared_ptr<const BoundaryCondition> >
NonlinearVariationalProblem::bcs() const
{
  return _bcs;
}
//-----------------------------------------------------------------------------
boost::shared_ptr<const FunctionSpace>
NonlinearVariationalProblem::trial_space() const
{
  dolfin_assert(_u);
  return _u->function_space();
}
//-----------------------------------------------------------------------------
boost::shared_ptr<const FunctionSpace>
NonlinearVariationalProblem::test_space() const
{
  dolfin_assert(_F);
  return _F->function_space(0);
}
//-----------------------------------------------------------------------------
bool NonlinearVariationalProblem::has_jacobian() const
{
  return _J; // cast to bool
}
//-----------------------------------------------------------------------------
void NonlinearVariationalProblem::check_forms() const
{
  // Check rank of residual F
  dolfin_assert(_F);
  if (_F->rank() != 1)
    dolfin_error("NonlinearVariationalProblem.cpp",
                 "define nonlinear variational problem F(u; v) = 0 for all v",
                 "Expecting the residual F to be a linear form (not rank %d)",
                 _F->rank());

  // Check rank of Jacobian J
  if (_J && _J->rank() != 2)
    dolfin_error("NonlinearVariationalProblem.cpp",
                 "define nonlinear variational problem F(u; v) = 0 for all v",
                 "Expecting the Jacobian J to be a bilinear form (not rank %d)",
                 _J->rank());

  /*
  // Check value of right-hand side
  if (rhs != 0)
    dolfin_error("NonlinearVariationalProblem.cpp",
                 "define nonlinear variational problem F(u; v) = 0 for all v",
                 "Expecting the right-hand side to be zero (not %d)",
                 rhs);
  */

  // FIXME: Should we add a check here that matches the function space
  // FIXME: of the solution variable u to a coefficient space for F?
  dolfin_assert(_u);
}
//-----------------------------------------------------------------------------
