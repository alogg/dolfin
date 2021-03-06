// Auto generated SWIG file for Python interface of DOLFIN
//
// Copyright (C) 2012 Kristian B. Oelgaard
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

// Autogenerated docstrings file, extracted from the DOLFIN source C++ files.

// Documentation extracted from: (module=adaptivity, header=GenericAdaptiveVariationalSolver.h)
%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver "
An abstract class for goal-oriented adaptive solution of
variational problems.

";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::solve "
Solve such that the functional error is less than the given
tolerance. Note that each call to solve is based on the
leaf-node of the variational problem

*Arguments*
    tol (float)
        The error tolerance
";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::solve_primal "
Solve the primal problem. Must be overloaded in subclass.

*Returns*
    :py:class:`Function`
        The solution to the primal problem
";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::extract_bcs "
Extract the boundary conditions for the primal problem. Must
be overloaded in subclass.

*Returns*
    list of :py:class:`BoundaryCondition`
        The primal boundary conditions
";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::evaluate_goal "
Evaluate the goal functional. Must be overloaded in subclass.

*Arguments*
   M (:py:class:`Form`)
       The functional to be evaluated
   u (:py:class:`Function`)
       The function of which to evaluate the functional

*Returns*
    float
        The value of M evaluated at u
";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::adapt_problem "
Adapt the problem to other mesh. Must be overloaded in subclass.

*Arguments*
   mesh (:py:class:`Mesh`)
       The other mesh
";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::adaptive_data "
Return stored adaptive data

*Returns*
   list of :py:class:`Parameters`
       The data stored in the adaptive loop
";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::default_parameters "
Default parameter values:

    \"max_iterations\"     (int)
    \"max_dimension\"      (int)
    \"plot_mesh\"          (bool)
    \"save_data\"          (bool)
    \"data_label\"         (std::string)
    \"reference\"          (double)
    \"marking_strategy\"   (std::string)
    \"marking_fraction\"   (double)
";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::summary "
Present summary of all adaptive data and parameters
";

%feature("docstring")  dolfin::GenericAdaptiveVariationalSolver::num_dofs_primal "
Return the number of degrees of freedom for primal problem

*Returns*
    _std::size_t_
        The number of degrees of freedom
";

// Documentation extracted from: (module=adaptivity, header=AdaptiveLinearVariationalSolver.h)
%feature("docstring")  dolfin::AdaptiveLinearVariationalSolver "
A class for goal-oriented adaptive solution of linear
variational problems.

For a linear variational problem of the form: find u in V
satisfying

    a(u, v) = L(v) for all v in :math:`\hat V`

and a corresponding conforming discrete problem: find u_h in V_h
satisfying

    a(u_h, v) = L(v) for all v in :math:`\hat V_h`

and a given goal functional M and tolerance tol, the aim is to
find a V_H and a u_H in V_H satisfying the discrete problem such
that

    \|M(u) - M(u_H)\| < tol

This strategy is based on dual-weighted residual error
estimators designed and automatically generated for the primal
problem and subsequent h-adaptivity.
";

%feature("docstring")  dolfin::AdaptiveLinearVariationalSolver::AdaptiveLinearVariationalSolver "
**Overloaded versions**

* AdaptiveLinearVariationalSolver\ (problem, goal)

  Create AdaptiveLinearVariationalSolver
  
  *Arguments*
      problem (:py:class:`LinearVariationalProblem`)
          The primal problem
      goal (:py:class:`GoalFunctional`)
          The goal functional

* AdaptiveLinearVariationalSolver\ (problem, goal)

  Create AdaptiveLinearVariationalSolver (shared ptr version)
  
  *Arguments*
      problem (:py:class:`LinearVariationalProblem`)
          The primal problem
      goal (:py:class:`GoalFunctional`)
          The goal functional

* AdaptiveLinearVariationalSolver\ (problem, goal, control)

  Create AdaptiveLinearVariationalSolver from variational
  problem, goal form and error control instance
  
  *Arguments*
      problem (:py:class:`LinearVariationalProblem`)
          The primal problem
      goal (:py:class:`Form`)
          The goal functional
      control (:py:class:`ErrorControl`)
          An error controller object
";

%feature("docstring")  dolfin::AdaptiveLinearVariationalSolver::solve_primal "
Solve the primal problem.

*Returns*
    :py:class:`Function`
        The solution to the primal problem
";

%feature("docstring")  dolfin::AdaptiveLinearVariationalSolver::extract_bcs "
Extract the boundary conditions for the primal problem.

*Returns*
    list of :py:class:`BoundaryCondition`
        The primal boundary conditions
";

%feature("docstring")  dolfin::AdaptiveLinearVariationalSolver::evaluate_goal "
Evaluate the goal functional.

*Arguments*
   M (:py:class:`Form`)
       The functional to be evaluated
   u (:py:class:`Function`)
       The function at which to evaluate the functional

*Returns*
    float
        The value of M evaluated at u
";

%feature("docstring")  dolfin::AdaptiveLinearVariationalSolver::adapt_problem "
Adapt the problem to other mesh.

*Arguments*
   mesh (:py:class:`Mesh`)
       The other mesh
";

%feature("docstring")  dolfin::AdaptiveLinearVariationalSolver::num_dofs_primal "
Return the number of degrees of freedom for primal problem

*Returns*
    _std::size_t_
        The number of degrees of freedom
";

%feature("docstring")  dolfin::AdaptiveLinearVariationalSolver::init "
Helper function for instance initialization

*Arguments*
   problem (:py:class:`LinearVariationalProblem`)
       The primal problem
   u (:py:class:`GoalFunctional`)
       The goal functional
";

// Documentation extracted from: (module=adaptivity, header=AdaptiveNonlinearVariationalSolver.h)
%feature("docstring")  dolfin::AdaptiveNonlinearVariationalSolver "
A class for goal-oriented adaptive solution of nonlinear
variational problems.

For a nonlinear variational problem of the form: find u in V
satisfying

    F(u; v) = 0 for all v in :math:`\hat V`

and a corresponding conforming discrete problem: find u_h in V_h
satisfying (at least approximately)

    F(u_h; v) = 0 for all v in :math:`\hat V_h`

and a given goal functional M and tolerance tol, the aim is to
find a V_H and a u_H in V_H satisfying the discrete problem such
that

    \|M(u) - M(u_H)\| < tol

This strategy is based on dual-weighted residual error
estimators designed and automatically generated for the primal
problem and subsequent h-adaptivity.
";

%feature("docstring")  dolfin::AdaptiveNonlinearVariationalSolver::AdaptiveNonlinearVariationalSolver "
**Overloaded versions**

* AdaptiveNonlinearVariationalSolver\ (problem, goal)

  Create AdaptiveNonlinearVariationalSolver
  
  *Arguments*
      problem (:py:class:`NonlinearVariationalProblem`)
          The primal problem
      goal (:py:class:`GoalFunctional`)
          The goal functional

* AdaptiveNonlinearVariationalSolver\ (problem, goal)

  Create AdaptiveNonlinearVariationalSolver (shared ptr version)
  
  *Arguments*
      problem (:py:class:`NonlinearVariationalProblem`)
          The primal problem
      goal (:py:class:`GoalFunctional`)
          The goal functional

* AdaptiveNonlinearVariationalSolver\ (problem, goal, control)

  Create AdaptiveLinearVariationalSolver from variational
  problem, goal form and error control instance
  
  *Arguments*
      problem (:py:class:`NonlinearVariationalProblem`)
          The primal problem
      goal (:py:class:`Form`)
          The goal functional
      control (:py:class:`ErrorControl`)
          An error controller object
";

%feature("docstring")  dolfin::AdaptiveNonlinearVariationalSolver::solve_primal "
Solve the primal problem.

*Returns*
    :py:class:`Function`
        The solution to the primal problem
";

%feature("docstring")  dolfin::AdaptiveNonlinearVariationalSolver::extract_bcs "
Extract the boundary conditions for the primal problem.

*Returns*
    list of :py:class:`BoundaryCondition`
        The primal boundary conditions
";

%feature("docstring")  dolfin::AdaptiveNonlinearVariationalSolver::evaluate_goal "
Evaluate the goal functional.

*Arguments*
   M (:py:class:`Form`)
       The functional to be evaluated
   u (:py:class:`Function`)
       The function at which to evaluate the functional

*Returns*
    float
        The value of M evaluated at u
";

%feature("docstring")  dolfin::AdaptiveNonlinearVariationalSolver::adapt_problem "
Adapt the problem to other mesh.

*Arguments*
   mesh (:py:class:`Mesh`)
       The other mesh
";

%feature("docstring")  dolfin::AdaptiveNonlinearVariationalSolver::num_dofs_primal "
Return the number of degrees of freedom for primal problem

*Returns*
    _std::size_t_
        The number of degrees of freedom
";

%feature("docstring")  dolfin::AdaptiveNonlinearVariationalSolver::init "
Helper function for instance initialization

*Arguments*
   problem (:py:class:`NonlinearVariationalProblem`)
       The primal problem
   u (:py:class:`GoalFunctional`)
       The goal functional
";

// Documentation extracted from: (module=adaptivity, header=GoalFunctional.h)
%feature("docstring")  dolfin::GoalFunctional "
A :py:class:`GoalFunctional` is a :py:class:`Form` of rank 0 with an associated
:py:class:`ErrorControl`.
";

%feature("docstring")  dolfin::GoalFunctional::GoalFunctional "
Create :py:class:`GoalFunctional`

*Arguments*
    rank (int)
        the rank of the functional (should be 0)
    num_coefficients (int)
        the number of coefficients in functional
";

%feature("docstring")  dolfin::GoalFunctional::update_ec "
Update error control instance with given forms

*Arguments*
    a (:py:class:`Form`)
        a bilinear form
    L (:py:class:`Form`)
        a linear form
";

// Documentation extracted from: (module=adaptivity, header=ErrorControl.h)
%feature("docstring")  dolfin::ErrorControl "
(Goal-oriented) Error Control class.
The notation used here follows the notation in \"Automated
goal-oriented error control I: stationary variational problems\",
ME Rognes and A Logg, 2010-2011.
";

%feature("docstring")  dolfin::ErrorControl::ErrorControl "
Create error control object

*Arguments*
    a_star (:py:class:`Form`)
       the bilinear form for the dual problem
    L_star (:py:class:`Form`)
       the linear form for the dual problem
    residual (:py:class:`Form`)
       a functional for the residual (error estimate)
    a_R_T (:py:class:`Form`)
       the bilinear form for the strong cell residual problem
    L_R_T (:py:class:`Form`)
       the linear form for the strong cell residual problem
    a_R_dT (:py:class:`Form`)
       the bilinear form for the strong facet residual problem
    L_R_dT (:py:class:`Form`)
       the linear form for the strong facet residual problem
    eta_T (:py:class:`Form`)
       a linear form over DG_0 for error indicators
    is_linear (bool)
       true iff primal problem is linear
";

%feature("docstring")  dolfin::ErrorControl::default_parameters "
Default parameter values:
";

%feature("docstring")  dolfin::ErrorControl::estimate_error "
Estimate the error relative to the goal M of the discrete
approximation 'u' relative to the variational formulation by
evaluating the weak residual at an approximation to the dual
solution.

*Arguments*
    u (:py:class:`Function`)
       the primal approximation

    bcs (list of :py:class:`BoundaryCondition`)
        the primal boundary conditions

*Returns*
    float
        error estimate
";

%feature("docstring")  dolfin::ErrorControl::compute_indicators "
Compute error indicators

*Arguments*
    indicators (:py:class:`MeshFunction`)
        the error indicators (to be computed)

    u (:py:class:`Function`)
        the primal approximation
";

%feature("docstring")  dolfin::ErrorControl::residual_representation "
Compute strong representation (strong cell and facet
residuals) of the weak residual.

*Arguments*
    R_T (:py:class:`Function`)
        the strong cell residual (to be computed)

    R_dT (:py:class:`SpecialFacetFunction`)
        the strong facet residual (to be computed)

    u (:py:class:`Function`)
        the primal approximation
";

%feature("docstring")  dolfin::ErrorControl::compute_cell_residual "
Compute representation for the strong cell residual
from the weak residual

*Arguments*
    R_T (:py:class:`Function`)
        the strong cell residual (to be computed)

    u (:py:class:`Function`)
        the primal approximation
";

%feature("docstring")  dolfin::ErrorControl::compute_facet_residual "
Compute representation for the strong facet residual from the
weak residual and the strong cell residual

*Arguments*
    R_dT (:py:class:`SpecialFacetFunction`)
        the strong facet residual (to be computed)

    u (:py:class:`Function`)
        the primal approximation

    R_T (:py:class:`Function`)
        the strong cell residual
";

%feature("docstring")  dolfin::ErrorControl::compute_dual "
Compute dual approximation defined by dual variational
problem and dual boundary conditions given by homogenized primal
boundary conditions.

*Arguments*
    z (:py:class:`Function`)
        the dual approximation (to be computed)

    bcs (list of :py:class:`BoundaryCondition`)
        the primal boundary conditions
";

%feature("docstring")  dolfin::ErrorControl::compute_extrapolation "
Compute extrapolation with boundary conditions

*Arguments*
    z (:py:class:`Function`)
        the extrapolated function (to be computed)

    bcs (list of :py:class:`BoundaryCondition`)
        the dual boundary conditions
";

// Documentation extracted from: (module=adaptivity, header=Extrapolation.h)
%feature("docstring")  dolfin::Extrapolation "
This class implements an algorithm for extrapolating a function
on a given function space from an approximation of that function
on a possibly lower-order function space.

This can be used to obtain a higher-order approximation of a
computed dual solution, which is necessary when the computed
dual approximation is in the test space of the primal problem,
thereby being orthogonal to the residual.

It is assumed that the extrapolation is computed on the same
mesh as the original function.
";

%feature("docstring")  dolfin::Extrapolation::extrapolate "
Compute extrapolation w from v
";

// Documentation extracted from: (module=adaptivity, header=LocalAssembler.h)
%feature("docstring")  dolfin::LocalAssembler "

";

%feature("docstring")  dolfin::LocalAssembler::assemble "

";

%feature("docstring")  dolfin::LocalAssembler::assemble_cell "

";

%feature("docstring")  dolfin::LocalAssembler::assemble_exterior_facet "

";

%feature("docstring")  dolfin::LocalAssembler::assemble_interior_facet "

";

// Documentation extracted from: (module=adaptivity, header=TimeSeries.h)
%feature("docstring")  dolfin::TimeSeries "
This class stores a time series of objects to file(s) in a
binary format which is efficient for reading and writing.

When objects are retrieved, the object stored at the time
closest to the given time will be used.

A new time series will check if values have been stored to
file before (for a series with the same name) and in that
case reuse those values. If new values are stored, old
values will be cleared.
";

%feature("docstring")  dolfin::TimeSeries::TimeSeries "
Create empty time series

*Arguments*
    name (str)
        The time series name
    compressed (bool)
        Use compressed file format (default false)
    store_connectivity (bool)
        Store all computed connectivity (default false)
";

%feature("docstring")  dolfin::TimeSeries::store "
**Overloaded versions**

* store\ (vector, t)

  Store vector at given time
  
  *Arguments*
      vector (:py:class:`GenericVector`)
          The vector to be stored.
      t (float)
          The time.

* store\ (mesh, t)

  Store mesh at given time
  
  *Arguments*
      mesh (:py:class:`Mesh`)
          The mesh to be stored.
      t (float)
          The time.
";

%feature("docstring")  dolfin::TimeSeries::retrieve "
**Overloaded versions**

* retrieve\ (vector, t, interpolate=true)

  Retrieve vector at given time
  
  *Arguments*
      vector (:py:class:`GenericVector`)
          The vector (values to be retrieved).
      t (float)
          The time.
      interpolate (bool)
          Optional argument: If true (default), interpolate
          time samples closest to t if t is not present.

* retrieve\ (mesh, t)

  Retrieve mesh at given time
  
  *Arguments*
      mesh (:py:class:`Mesh`)
          The mesh (values to be retrieved).
      t (float)
          The time.
";

%feature("docstring")  dolfin::TimeSeries::vector_times "
Return array of sample times for vectors

*Returns*
    numpy.array(float)
        The times.
";

%feature("docstring")  dolfin::TimeSeries::mesh_times "
Return array of sample times for meshes

*Returns*
    numpy.array(float)
        The times.
";

%feature("docstring")  dolfin::TimeSeries::clear "
Clear time series
";

%feature("docstring")  dolfin::TimeSeries::filename_data "
Return filename for data

*Arguments*
    series_name (str)
        The time series name
    type_name (str)
        The type of data
    index (std::size_t)
        The index
    compressed (bool)
        True if compressed file format

*Returns*
    str
        The filename
";

%feature("docstring")  dolfin::TimeSeries::filename_times "
Return filename for times

*Arguments*
    series_name (str)
        The time series name
    type_name (str)
        The type of data
    compressed (bool)
        True if compressed file format

*Returns*
    str
        The filename
";

%feature("docstring")  dolfin::TimeSeries::str "
Return informal string representation (pretty-print)
";

%feature("docstring")  dolfin::TimeSeries::default_parameters "
Default parameter values
";

// Documentation extracted from: (module=adaptivity, header=adapt.h)
%feature("docstring")  dolfin::adapt "
**Overloaded versions**

* adapt\ (mesh)

  Refine mesh uniformly

* adapt\ (mesh, cell_markers)

  Refine mesh based on cell markers

* adapt\ (space)

  Refine function space uniformly

* adapt\ (space, cell_markers)

  Refine function space based on cell markers

* adapt\ (space, adapted_mesh)

  Refine function space based on refined mesh

* adapt\ (function, adapted_mesh, interpolate=true)

  Adapt Function based on adapted mesh
  
  *Arguments*
      function (:py:class:`Function`)
          The function that should be adapted
      adapted_mesh (:py:class:`Mesh`)
          The new mesh
      interpolate (bool)
          Optional argument, default is true. If false, the
          function's function space is adapted, but the values are
          not interpolated.
  
  *Returns*
      :py:class:`Function`
          The adapted function

* adapt\ (function, adapted_mesh)

  Refine GenericFunction based on refined mesh

* adapt\ (mesh_function, adapted_mesh)

  Refine mesh function<std::size_t> based on mesh

* adapt\ (bc, adapted_mesh, S)

  Refine Dirichlet bc based on refined mesh

* adapt\ (form, adapted_mesh, adapt_coefficients=true)

  Adapt form based on adapted mesh
  
  *Arguments*
      form (:py:class:`Form`)
          The form that should be adapted
      adapted_mesh (:py:class:`Mesh`)
          The new mesh
      adapt_coefficients (bool)
          Optional argument, default is true. If false, the form
          coefficients are not explictly adapted, but pre-adapted
          coefficients will be transferred.
  
  *Returns*
      :py:class:`Form`
          The adapted form

* adapt\ (problem, adapted_mesh)

  Refine linear variational problem based on mesh

* adapt\ (problem, adapted_mesh)

  Refine nonlinear variational problem based on mesh

* adapt\ (ec, adapted_mesh, adapt_coefficients=true)

  Adapt error control object based on adapted mesh
  
  *Arguments*
      ec (:py:class:`ErrorControl`)
          The error control object to be adapted
      adapted_mesh (:py:class:`Mesh`)
          The new mesh
      adapt_coefficients (bool)
          Optional argument, default is true. If false, any form
          coefficients are not explictly adapted, but pre-adapted
          coefficients will be transferred.
  
  *Returns*
      :py:class:`ErrorControl`
          The adapted error control object
";

%feature("docstring")  dolfin::adapt_markers "
Helper function for refinement of boundary conditions
";

// Documentation extracted from: (module=adaptivity, header=marking.h)
%feature("docstring")  dolfin::mark "
Mark cells based on indicators and given marking strategy

*Arguments*
    markers (:py:class:`MeshFunction`)
        the cell markers (to be computed)

    indicators (:py:class:`MeshFunction`)
        error indicators (one per cell)

    strategy (str)
        the marking strategy

    fraction (float)
        the marking fraction
";

%feature("docstring")  dolfin::dorfler_mark "
Mark cells using Dorfler marking

*Arguments*
    markers (:py:class:`MeshFunction`)
        the cell markers (to be computed)

    indicators (:py:class:`MeshFunction`)
        error indicators (one per cell)

    fraction (float)
        the marking fraction
";

// Documentation extracted from: (module=adaptivity, header=adaptivesolve.h)
%feature("docstring")  dolfin::solve "
**Overloaded versions**

* solve\ (equation, u, tol, M)

  Solve linear variational problem a(u, v) == L(v) without
  essential boundary conditions

* solve\ (equation, u, bc, tol, M)

  Solve linear variational problem a(u, v) == L(v) with single
  boundary condition

* solve\ (equation, u, bcs, tol, M)

  Solve linear variational problem a(u, v) == L(v) with list of
  boundary conditions

* solve\ (equation, u, J, tol, M)

  Solve nonlinear variational problem F(u; v) = 0 without
  essential boundary conditions

* solve\ (equation, u, bc, J, tol, M)

  Solve linear variational problem F(u; v) = 0 with single
  boundary condition

* solve\ (equation, u, bcs, J, tol, M)

  Solve linear variational problem F(u; v) = 0 with list of
  boundary conditions
";

