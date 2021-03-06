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

// Documentation extracted from: (module=function, header=GenericFunction.h)
%feature("docstring")  dolfin::GenericFunction "
This is a common base class for functions. Functions can be
evaluated at a given point and they can be restricted to a given
cell in a finite element mesh. This functionality is implemented
by sub-classes that implement the eval() and restrict() functions.

DOLFIN provides two implementations of the GenericFunction
interface in the form of the classes Function and Expression.

Sub-classes may optionally implement the update() function that
will be called prior to restriction when running in parallel.
";

%feature("docstring")  dolfin::GenericFunction::GenericFunction "
Constructor
";

%feature("docstring")  dolfin::GenericFunction::value_rank "
Return value rank
";

%feature("docstring")  dolfin::GenericFunction::value_dimension "
Return value dimension for given axis
";

%feature("docstring")  dolfin::GenericFunction::eval "
**Overloaded versions**

* eval\ (values, x, cell)

  Evaluate at given point in given cell

* eval\ (values, x)

  Evaluate at given point
";

%feature("docstring")  dolfin::GenericFunction::restrict "
Restrict function to local cell (compute expansion coefficients w)
";

%feature("docstring")  dolfin::GenericFunction::compute_vertex_values "
Compute values at all mesh vertices
";

%feature("docstring")  dolfin::GenericFunction::update "
Update off-process ghost coefficients
";

%feature("docstring")  dolfin::GenericFunction::operator "
**Overloaded versions**

* operator\ (x)

  Evaluation at given point (scalar function)

* operator\ (x, y)

  Evaluation at given point (scalar function)

* operator\ (x, y, z)

  Evaluation at given point (scalar function)

* operator\ (p)

  Evaluation at given point (scalar function)

* operator\ (values, x)

  Evaluation at given point (vector-valued function)

* operator\ (values, x, y)

  Evaluation at given point (vector-valued function)

* operator\ (values, x, y, z)

  Evaluation at given point (vector-valued function)

* operator\ (values, p)

  Evaluation at given point (vector-valued function)
";

%feature("docstring")  dolfin::GenericFunction::value_size "
Evaluation at given point
Return value size (product of value dimensions)
";

%feature("docstring")  dolfin::GenericFunction::evaluate "
Evaluate function at given point in cell
";

// Documentation extracted from: (module=function, header=Expression.h)
%feature("docstring")  dolfin::Expression "
This class represents a user-defined expression. Expressions can
be used as coefficients in variational forms or interpolated
into finite element spaces.

An expression is defined by overloading the eval() method. Users
may choose to overload either a simple version of eval(), in the
case of expressions only depending on the coordinate x, or an
optional version for expressions depending on x and mesh data
like cell indices or facet normals.

The geometric dimension (the size of x) and the value rank and
dimensions of an expression must supplied as arguments to the
constructor.
";

%feature("docstring")  dolfin::Expression::Expression "
**Overloaded versions**

* Expression\ ()

  Create scalar expression.

* Expression\ (dim)

  Create vector-valued expression with given dimension.
  
  *Arguments*
      dim (std::size_t)
          Dimension of the vector-valued expression.

* Expression\ (dim0, dim1)

  Create matrix-valued expression with given dimensions.
  
  *Arguments*
      dim0 (std::size_t)
          Dimension (rows).
      dim1 (std::size_t)
          Dimension (columns).

* Expression\ (value_shape)

  Create tensor-valued expression with given shape.
  
  *Arguments*
      value_shape (std::vector<std::size_t>)
          Shape of expression.

* Expression\ (expression)

  Copy constructor
  
  *Arguments*
      expression (:py:class:`Expression`)
          Object to be copied.
";

%feature("docstring")  dolfin::Expression::eval "
**Overloaded versions**

* eval\ (values, x, cell)

  Note: The reimplementation of eval is needed for the Python interface.
  Evaluate at given point in given cell.
  
  *Arguments*
      values (numpy.array(float))
          The values at the point.
      x (numpy.array(float))
          The coordinates of the point.
      cell (ufc::cell)
          The cell which contains the given point.

* eval\ (values, x)

  Evaluate at given point.
  
  *Arguments*
      values (numpy.array(float))
          The values at the point.
      x (numpy.array(float))
          The coordinates of the point.
";

%feature("docstring")  dolfin::Expression::value_rank "
Return value rank.

*Returns*
    std::size_t
        The value rank.
";

%feature("docstring")  dolfin::Expression::value_dimension "
Return value dimension for given axis.

*Arguments*
    i (std::size_t)
        Integer denoting the axis to use.

*Returns*
    std::size_t
        The value dimension (for the given axis).
";

%feature("docstring")  dolfin::Expression::restrict "
Restrict function to local cell (compute expansion coefficients w).

*Arguments*
    w (list of doubles)
        Expansion coefficients.
    element (:py:class:`FiniteElement`)
        The element.
    dolfin_cell (:py:class:`Cell`)
        The cell.
    ufc_cell (ufc::cell)
        The ufc::cell.
";

%feature("docstring")  dolfin::Expression::compute_vertex_values "
Compute values at all mesh vertices.

*Arguments*
    vertex_values (numpy.array(float))
        The values at all vertices.
    mesh (:py:class:`Mesh`)
        The mesh.
";

// Documentation extracted from: (module=function, header=Function.h)
%feature("docstring")  dolfin::Function "
This class represents a function :math:`u_h` in a finite
element function space :math:`V_h`, given by

.. math::

    u_h = \sum_{i=1}^{n} U_i \phi_i

where :math:`\{\phi_i\}_{i=1}^{n}` is a basis for :math:`V_h`,
and :math:`U` is a vector of expansion coefficients for :math:`u_h`.
";

%feature("docstring")  dolfin::Function::Function "
**Overloaded versions**

* Function\ (V)

  Create function on given function space
  
  *Arguments*
      V (:py:class:`FunctionSpace`)
          The function space.
  
  *Example*
      .. note::
      
          No example code available for this function.

* Function\ (V)

  Create function on given function space (shared data)
  
  *Arguments*
      V (:py:class:`FunctionSpace`)
          The function space.

* Function\ (V, x)

  Create function on given function space with a given vector
  (shared data)
  
  *Warning: This constructor is intended for internal library use only*
  
  *Arguments*
      V (:py:class:`FunctionSpace`)
          The function space.
      x (:py:class:`GenericVector`)
          The vector.

* Function\ (V, filename)

  Create function from vector of dofs stored to file
  
  *Arguments*
      V (:py:class:`FunctionSpace`)
          The function space.
      filename_vector (str)
          The name of the file containing the vector.
      filename_dofdata (str)
          The name of the file containing the dofmap data.

* Function\ (V, filename)

  Create function from vector of dofs stored to file (shared data)
  
  *Arguments*
      V (:py:class:`FunctionSpace`)
          The function space.
      filename_dofdata (str)
          The name of the file containing the dofmap data.

* Function\ (v)

  Copy constructor
  
  *Arguments*
      v (:py:class:`Function`)
          The object to be copied.

* Function\ (v, i)

  Sub-function constructor with shallow copy of vector (used in Python
  interface)
  
  *Arguments*
      v (:py:class:`Function`)
          The function to be copied.
      i (std::size_t)
          Index of subfunction.
  
";

%feature("docstring")  dolfin::Function::operator= "
**Overloaded versions**

* operator=\ (v)

  Assignment from function
  
  *Arguments*
      v (:py:class:`Function`)
          Another function.

* operator=\ (v)

  Assignment from expression using interpolation
  
  *Arguments*
      v (:py:class:`Expression`)
          The expression.
";

%feature("docstring")  dolfin::Function::operator[] "
Extract subfunction

*Arguments*
    i (std::size_t)
        Index of subfunction.
";

%feature("docstring")  dolfin::Function::function_space "
Return shared pointer to function space

*Returns*
    :py:class:`FunctionSpace`
        Return the shared pointer.
";

%feature("docstring")  dolfin::Function::vector "
**Overloaded versions**

* vector\ ()

  Return vector of expansion coefficients (non-const version)
  
  *Returns*
      :py:class:`GenericVector`
          The vector of expansion coefficients.

* vector\ ()

  Return vector of expansion coefficients (const version)
  
  *Returns*
      :py:class:`GenericVector`
          The vector of expansion coefficients (const).
";

%feature("docstring")  dolfin::Function::in "
Check if function is a member of the given function space

*Arguments*
    V (:py:class:`FunctionSpace`)
        The function space.

*Returns*
    bool
        True if the function is in the function space.
";

%feature("docstring")  dolfin::Function::geometric_dimension "
Return geometric dimension

*Returns*
    std::size_t
        The geometric dimension.
";

%feature("docstring")  dolfin::Function::eval "
**Overloaded versions**

* eval\ (values, x)

  Evaluate function at given coordinates
  
  *Arguments*
      values (numpy.array(float))
          The values.
      x (numpy.array(float))
          The coordinates.

* eval\ (values, x, dolfin_cell, ufc_cell)

  Evaluate function at given coordinates in given cell
  
  *Arguments*
      values (numpy.array(float))
          The values.
      x (numpy.array(float))
          The coordinates.
      dolfin_cell (:py:class:`Cell`)
          The cell.
      ufc_cell (ufc::cell)
          The ufc::cell.

* eval\ (values, x, cell)

  Evaluate at given point in given cell
  
  *Arguments*
      values (numpy.array(float))
          The values at the point.
      x (numpy.array(float))
          The coordinates of the point.
      cell (ufc::cell)
          The cell which contains the given point.
";

%feature("docstring")  dolfin::Function::interpolate "
Interpolate function (on possibly non-matching meshes)

*Arguments*
    v (:py:class:`GenericFunction`)
        The function to be interpolated.
";

%feature("docstring")  dolfin::Function::extrapolate "
Extrapolate function (from a possibly lower-degree function space)

*Arguments*
    v (:py:class:`Function`)
        The function to be extrapolated.
";

%feature("docstring")  dolfin::Function::value_rank "
Return value rank

*Returns*
    std::size_t
        The value rank.
";

%feature("docstring")  dolfin::Function::value_dimension "
Return value dimension for given axis

*Arguments*
    i (std::size_t)
        The index of the axis.

*Returns*
    std::size_t
        The value dimension.
";

%feature("docstring")  dolfin::Function::non_matching_eval "
Evaluate function for given data (non-matching meshes)

*Arguments*
    values (numpy.array(float))
        The values at the point.
    x (numpy.array(float))
        The coordinates of the point.
    cell (ufc::cell)
        The cell.
";

%feature("docstring")  dolfin::Function::restrict "
Restrict function to local cell (compute expansion coefficients w)

*Arguments*
    w (list of doubles)
        Expansion coefficients.
    element (:py:class:`FiniteElement`)
        The element.
    dolfin_cell (:py:class:`Cell`)
        The cell.
    ufc_cell (ufc::cell)
        The ufc::cell.
";

%feature("docstring")  dolfin::Function::compute_vertex_values "
**Overloaded versions**

* compute_vertex_values\ (vertex_values, mesh)

  Compute values at all mesh vertices
  
  *Arguments*
      vertex_values (numpy.array(float))
          The values at all vertices.
      mesh (:py:class:`Mesh`)
          The mesh.

* compute_vertex_values\ (vertex_values)

  Compute values at all mesh vertices
  
  *Arguments*
      vertex_values (numpy.array(float))
          The values at all vertices.
";

%feature("docstring")  dolfin::Function::update "
Update off-process ghost coefficients
";

// Documentation extracted from: (module=function, header=FunctionSpace.h)
%feature("docstring")  dolfin::FunctionSpace "
This class represents a finite element function space defined by
a mesh, a finite element, and a local-to-global mapping of the
degrees of freedom (dofmap).
";

%feature("docstring")  dolfin::FunctionSpace::FunctionSpace "
**Overloaded versions**

* FunctionSpace\ (mesh, element, dofmap)

  Create function space for given mesh, element and dofmap
  (shared data)
  
  *Arguments*
      mesh (:py:class:`Mesh`)
          The mesh.
      element (:py:class:`FiniteElement`)
          The element.
      dofmap (:py:class:`GenericDofMap`)
          The dofmap.

* FunctionSpace\ (mesh)

  Create empty function space for later initialization. This
  constructor is intended for use by any sub-classes which need
  to construct objects before the initialisation of the base
  class. Data can be attached to the base class using
  FunctionSpace::attach(...).
  
  *Arguments*
      mesh (:py:class:`Mesh`)
          The mesh.

* FunctionSpace\ (V)

  Copy constructor
  
  *Arguments*
      V (:py:class:`FunctionSpace`)
          The object to be copied.
";

%feature("docstring")  dolfin::FunctionSpace::attach "
Attach data to an empty function space

*Arguments*
    element (:py:class:`FiniteElement`)
        The element.
    dofmap (:py:class:`GenericDofMap`)
        The dofmap.
";

%feature("docstring")  dolfin::FunctionSpace::operator= "
Assignment operator

*Arguments*
    V (:py:class:`FunctionSpace`)
        Another function space.
";

%feature("docstring")  dolfin::FunctionSpace::operator== "
Equality operator

*Arguments*
    V (:py:class:`FunctionSpace`)
        Another function space.
";

%feature("docstring")  dolfin::FunctionSpace::operator!= "
Unequality operator

*Arguments*
    V (:py:class:`FunctionSpace`)
        Another function space.
";

%feature("docstring")  dolfin::FunctionSpace::mesh "
Return mesh

*Returns*
    :py:class:`Mesh`
        The mesh.
";

%feature("docstring")  dolfin::FunctionSpace::element "
Return finite element

*Returns*
    :py:class:`FiniteElement`
        The finite element.
";

%feature("docstring")  dolfin::FunctionSpace::dofmap "
Return dofmap

*Returns*
    :py:class:`GenericDofMap`
        The dofmap.
";

%feature("docstring")  dolfin::FunctionSpace::dim "
Return dimension of function space

*Returns*
    std::size_t
        The dimension of the function space.
";

%feature("docstring")  dolfin::FunctionSpace::interpolate "
Interpolate function v into function space, returning the
vector of expansion coefficients

*Arguments*
    expansion_coefficients (:py:class:`GenericVector`)
        The expansion coefficients.
    v (:py:class:`GenericFunction`)
        The function to be interpolated.
";

%feature("docstring")  dolfin::FunctionSpace::operator[] "
Extract subspace for component

*Arguments*
    i (std::size_t)
        Index of the subspace.
*Returns*
    :py:class:`FunctionSpace`
        The subspace.
";

%feature("docstring")  dolfin::FunctionSpace::extract_sub_space "
Extract subspace for component

*Arguments*
    component (std::vector<std::size_t>)
        The component.

*Returns*
    :py:class:`FunctionSpace`
        The subspace.
";

%feature("docstring")  dolfin::FunctionSpace::collapse "
**Overloaded versions**

* collapse\ ()

  Collapse a subspace and return a new function space
  
  *Returns*
      :py:class:`FunctionSpace`
          The new function space.

* collapse\ (collapsed_dofs)

  Collapse a subspace and return a new function space and a map
  from new to old dofs
  
  *Arguments*
      collapsed_dofs (boost::unordered_map<std::size_t, std::size_t>)
          The map from new to old dofs.
  
  *Returns*
      :py:class:`FunctionSpace`
        The new function space.
";

%feature("docstring")  dolfin::FunctionSpace::has_cell "
Check if function space has given cell

*Arguments*
    cell (:py:class:`Cell`)
        The cell.

*Returns*
    bool
        True if the function space has the given cell.
";

%feature("docstring")  dolfin::FunctionSpace::has_element "
Check if function space has given element

*Arguments*
    element (:py:class:`FiniteElement`)
        The finite element.

*Returns*
    bool
        True if the function space has the given element.
";

%feature("docstring")  dolfin::FunctionSpace::component "
Return component

*Returns*
    std::vector<std::size_t>
        The component (relative to superspace).
";

%feature("docstring")  dolfin::FunctionSpace::str "
Return informal string representation (pretty-print)

*Arguments*
    verbose (bool)
        Flag to turn on additional output.

*Returns*
    str
        An informal representation of the function space.
";

%feature("docstring")  dolfin::FunctionSpace::print_dofmap "
Print dofmap (useful for debugging)
";

// Documentation extracted from: (module=function, header=SubSpace.h)
%feature("docstring")  dolfin::SubSpace "
This class represents a subspace (component) of a function space.

The subspace is specified by an array of indices. For example,
the array [3, 0, 2] specifies subspace 2 of subspace 0 of
subspace 3.

A typical example is the function space W = V x P for Stokes.
Here, V = W[0] is the subspace for the velocity component and
P = W[1] is the subspace for the pressure component. Furthermore,
W[0][0] = V[0] is the first component of the velocity space etc.
";

%feature("docstring")  dolfin::SubSpace::SubSpace "
**Overloaded versions**

* SubSpace\ (V, component)

  Create subspace for given component (one level)

* SubSpace\ (V, component, sub_component)

  Create subspace for given component (two levels)

* SubSpace\ (V, component)

  Create subspace for given component (n levels)
";

// Documentation extracted from: (module=function, header=Constant.h)
%feature("docstring")  dolfin::Constant "
This class represents a constant-valued expression.
";

%feature("docstring")  dolfin::Constant::Constant "
**Overloaded versions**

* Constant\ (value)

  Create scalar constant
  
  *Arguments*
      value (float)
          The scalar to create a Constant object from.
  
  *Example*
      .. note::
      
          No example code available for this function.

* Constant\ (value0, value1)

  Create vector constant (dim = 2)
  
  *Arguments*
      value0 (float)
          The first vector element.
      value1 (float)
          The second vector element.
  
  *Example*
      .. note::
      
          No example code available for this function.

* Constant\ (value0, value1, value2)

  Create vector constant (dim = 3)
  
  *Arguments*
      value0 (float)
          The first vector element.
      value1 (float)
          The second vector element.
      value2 (float)
          The third vector element.
  
  *Example*
      .. note::
      
          No example code available for this function.

* Constant\ (values)

  Create vector-valued constant
  
  *Arguments*
      values (numpy.array(float))
          Values to create a vector-valued constant from.

* Constant\ (value_shape, values)

  Create tensor-valued constant for flattened array of values
  
  *Arguments*
      value_shape (std::vector<std::size_t>)
          Shape of tensor.
      values (numpy.array(float))
          Values to create tensor-valued constant from.

* Constant\ (constant)

  Copy constructor
  
  *Arguments*
      constant (:py:class:`Constant`)
          Object to be copied.
";

%feature("docstring")  dolfin::Constant::operator= "
**Overloaded versions**

* operator=\ (constant)

  Assignment operator
  
  *Arguments*
      constant (:py:class:`Constant`)
          Another constant.

* operator=\ (constant)

  Assignment operator
  
  *Arguments*
      constant (float)
          Another constant.
";

%feature("docstring")  dolfin::Constant::operator double "
Cast to double (for scalar constants)

*Returns*
    float
        The scalar value.
";

// Documentation extracted from: (module=function, header=SpecialFunctions.h)
%feature("docstring")  dolfin::MeshCoordinates "
This Function represents the mesh coordinates on a given mesh.
";

%feature("docstring")  dolfin::MeshCoordinates::MeshCoordinates "
Constructor
";

%feature("docstring")  dolfin::MeshCoordinates::eval "
Evaluate function
";

%feature("docstring")  dolfin::FacetArea "
This function represents the area/length of a cell facet on a given mesh.
";

%feature("docstring")  dolfin::FacetArea::FacetArea "
Constructor
";

%feature("docstring")  dolfin::FacetArea::eval "
Evaluate function
";

// Documentation extracted from: (module=function, header=SpecialFacetFunction.h)
%feature("docstring")  dolfin::SpecialFacetFunction::SpecialFacetFunction "
**Overloaded versions**

* SpecialFacetFunction\ (f_e)

  Create (scalar-valued) SpecialFacetFunction
  
  *Arguments*
      f_e (list of :py:class:`Function`)
         Separate _Function_s for each facet

* SpecialFacetFunction\ (f_e, dim)

  Create (vector-valued) SpecialFacetFunction
  
  *Arguments*
      f_e (list of :py:class:`Function`)
         Separate _Function_s for each facet
  
      dim (int)
          The value-dimension of the Functions
";

%feature("docstring")  dolfin::SpecialFacetFunction::eval "
Evaluate SpecialFacetFunction (cf :py:class:`Expression`.eval)
Evaluate function for given cell
";

%feature("docstring")  dolfin::SpecialFacetFunction::operator[] "
Extract sub-function i

*Arguments*
    i (int)
       component

*Returns*
    :py:class:`Function`
";

