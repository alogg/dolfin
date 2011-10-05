/* -*- C -*- */
// Copyright (C) 2007-2009 Johan Hake
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
// First added:  2008-11-02
// Last changed: 2011-10-05

//-----------------------------------------------------------------------------
// Extend FunctionSpace so one can check if a Function is in a FunctionSpace
//-----------------------------------------------------------------------------
%extend dolfin::FunctionSpace {
%pythoncode %{
def __contains__(self,u):
    "Check whether a function is in the FunctionSpace"
    assert(isinstance(u,Function))
    return u._in(self)
%}
}

//-----------------------------------------------------------------------------
// Extend Function so f.function_space() returns a dolfin.FunctionSpace
//-----------------------------------------------------------------------------
%extend dolfin::Function {
%pythoncode %{
def function_space(self):
    "Return the FunctionSpace"
    from dolfin.functions.functionspace import FunctionSpaceFromCpp
    return FunctionSpaceFromCpp(self._function_space())
%}
}


//-----------------------------------------------------------------------------
// Copy
//-----------------------------------------------------------------------------
%extend dolfin::Function {
%pythoncode %{
def copy(self, deepcopy=False):
    "Return a dolfin.Function of itself"
    from dolfin.functions.function import Function
    if deepcopy:
        return Function(self.function_space(), self.vector().copy())
    return Function(self.function_space(), self.vector())
%}
}

//-----------------------------------------------------------------------------
// Extend Function so f.leaf_node()/root_node() returns a dolfin.Function.
// Not doing this on purpose for child()/parent().
// -----------------------------------------------------------------------------
%extend dolfin::Function {
%pythoncode %{
def leaf_node(self):
    "Return the finest Function in hierarchy"
    f = HierarchicalFunction.leaf_node(self)
    return f.copy()
%}
}
%extend dolfin::Function {
%pythoncode %{
def root_node(self):
    "Return the coarsest Function in hierarchy"
    f = HierarchicalFunction.root_node(self)
    return f.copy()
%}
}
