// Auto generated SWIG file for Python interface of DOLFIN
//
// Copyright (C) 2012 Johan Hake
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


// The PyDOLFIN extension module for the mesh module
%module(package="dolfin.cpp.mesh", directors="1") mesh

// Define module name for conditional includes
#define MESHMODULE

%{

// Include types from dependent modules

// #include types from common submodule of module common
#include "dolfin/common/Array.h"
#include "dolfin/common/Variable.h"
#include "dolfin/common/Hierarchical.h"

// #include types from parameter submodule of module common
#include "dolfin/parameter/Parameter.h"
#include "dolfin/parameter/Parameters.h"

// #include types from la submodule of module la
#include "dolfin/la/LinearAlgebraObject.h"

// #include types from function submodule of module function
#include "dolfin/function/GenericFunction.h"
#include "dolfin/function/Expression.h"
#include "dolfin/function/Function.h"
#include "dolfin/function/FunctionSpace.h"

// #include types from quadrature submodule of module fem
#include "dolfin/quadrature/BarycenterQuadrature.h"

// #include types from fem submodule of module fem
#include "dolfin/fem/GenericDofMap.h"
#include "dolfin/fem/DofMap.h"

// Include types from present module mesh

// #include types from intersection submodule
#include "dolfin/intersection/IntersectionOperator.h"
#include "dolfin/intersection/PrimitiveIntersector.h"
#include "dolfin/intersection/PrimitiveTraits.h"
#include "dolfin/intersection/MeshPrimitive.h"

// #include types from mesh submodule
#include "dolfin/mesh/CellType.h"
#include "dolfin/mesh/MeshTopology.h"
#include "dolfin/mesh/MeshGeometry.h"
#include "dolfin/mesh/MeshDomains.h"
#include "dolfin/mesh/MeshData.h"
#include "dolfin/mesh/Mesh.h"
#include "dolfin/mesh/MeshEntity.h"
#include "dolfin/mesh/MeshEntityIterator.h"
#include "dolfin/mesh/MeshEntityIteratorBase.h"
#include "dolfin/mesh/SubsetIterator.h"
#include "dolfin/mesh/Point.h"
#include "dolfin/mesh/Vertex.h"
#include "dolfin/mesh/Edge.h"
#include "dolfin/mesh/Face.h"
#include "dolfin/mesh/Facet.h"
#include "dolfin/mesh/Cell.h"
#include "dolfin/mesh/FacetCell.h"
#include "dolfin/mesh/MeshConnectivity.h"
#include "dolfin/mesh/MeshEditor.h"
#include "dolfin/mesh/DynamicMeshEditor.h"
#include "dolfin/mesh/LocalMeshValueCollection.h"
#include "dolfin/mesh/MeshFunction.h"
#include "dolfin/mesh/MeshPartitioning.h"
#include "dolfin/mesh/MeshValueCollection.h"
#include "dolfin/mesh/MeshColoring.h"
#include "dolfin/mesh/MeshRenumbering.h"
#include "dolfin/mesh/MeshTransformation.h"
#include "dolfin/mesh/LocalMeshData.h"
#include "dolfin/mesh/SubDomain.h"
#include "dolfin/mesh/SubMesh.h"
#include "dolfin/mesh/Restriction.h"
#include "dolfin/mesh/DomainBoundary.h"
#include "dolfin/mesh/BoundaryMesh.h"
#include "dolfin/mesh/PeriodicBoundaryComputation.h"

// #include types from generation submodule
#include "dolfin/generation/PolygonalMeshGenerator.h"
#include "dolfin/generation/PolyhedralMeshGenerator.h"
#include "dolfin/generation/Triangulate.h"
#include "dolfin/generation/BoxMesh.h"
#include "dolfin/generation/IntervalMesh.h"
#include "dolfin/generation/Interval.h"
#include "dolfin/generation/RectangleMesh.h"
#include "dolfin/generation/UnitTetrahedronMesh.h"
#include "dolfin/generation/UnitCubeMesh.h"
#include "dolfin/generation/UnitCube.h"
#include "dolfin/generation/UnitIntervalMesh.h"
#include "dolfin/generation/UnitInterval.h"
#include "dolfin/generation/UnitTriangleMesh.h"
#include "dolfin/generation/UnitSquareMesh.h"
#include "dolfin/generation/UnitSquare.h"
#include "dolfin/generation/UnitCircleMesh.h"
#include "dolfin/generation/UnitCircle.h"
#include "dolfin/generation/CSGGeometry.h"
#include "dolfin/generation/CSGMeshGenerator.h"
#include "dolfin/generation/CSGCGALMeshGenerator2D.h"
#include "dolfin/generation/CSGCGALMeshGenerator3D.h"
#include "dolfin/generation/CSGOperators.h"
#include "dolfin/generation/CSGPrimitive.h"
#include "dolfin/generation/CSGPrimitives2D.h"
#include "dolfin/generation/CSGPrimitives3D.h"
#include "dolfin/generation/CSGGeometries3D.h"

// #include types from refinement submodule
#include "dolfin/refinement/refine.h"

// #include types from graph submodule
#include "dolfin/graph/Graph.h"
#include "dolfin/graph/GraphBuilder.h"
#include "dolfin/graph/BoostGraphOrdering.h"
#include "dolfin/graph/SCOTCH.h"

// #include types from ale submodule
#include "dolfin/ale/ALE.h"
#include "dolfin/ale/MeshDisplacement.h"

// NumPy includes
#define PY_ARRAY_UNIQUE_SYMBOL PyDOLFIN_MESH
#include <numpy/arrayobject.h>
%}

%init%{
import_array();
%}

// Include global SWIG interface files:
// Typemaps, shared_ptr declarations, exceptions, version
%include "dolfin/swig/globalincludes.i"

// %import types from submodule common of SWIG module common
%include "dolfin/swig/common/pre.i"
%import(module="common") "dolfin/common/Array.h"
%import(module="common") "dolfin/common/Variable.h"
%import(module="common") "dolfin/common/Hierarchical.h"

// %import types from submodule parameter of SWIG module common
%include "dolfin/swig/parameter/pre.i"
%import(module="common") "dolfin/parameter/Parameter.h"
%import(module="common") "dolfin/parameter/Parameters.h"

// %import types from submodule la of SWIG module la
%include "dolfin/swig/la/pre.i"
%import(module="la") "dolfin/la/LinearAlgebraObject.h"

// %import types from submodule function of SWIG module function
%include "dolfin/swig/function/pre.i"
%import(module="function") "dolfin/function/GenericFunction.h"
%import(module="function") "dolfin/function/Expression.h"
%import(module="function") "dolfin/function/Function.h"
%import(module="function") "dolfin/function/FunctionSpace.h"

// %import types from submodule quadrature of SWIG module fem
%import(module="fem") "dolfin/quadrature/BarycenterQuadrature.h"

// %import types from submodule fem of SWIG module fem
%include "dolfin/swig/fem/pre.i"
%import(module="fem") "dolfin/fem/GenericDofMap.h"
%import(module="fem") "dolfin/fem/DofMap.h"

// Turn on SWIG generated signature documentation and include doxygen
// generated docstrings
//%feature("autodoc", "1");
%include "dolfin/swig/intersection/docstrings.i"
%include "dolfin/swig/mesh/docstrings.i"
%include "dolfin/swig/generation/docstrings.i"
%include "dolfin/swig/refinement/docstrings.i"
%include "dolfin/swig/graph/docstrings.i"
%include "dolfin/swig/ale/docstrings.i"

// %include types from submodule intersection
%include "dolfin/intersection/IntersectionOperator.h"
%include "dolfin/intersection/PrimitiveIntersector.h"
%include "dolfin/intersection/PrimitiveTraits.h"
%include "dolfin/intersection/MeshPrimitive.h"

// %include types from submodule mesh
%include "dolfin/swig/mesh/pre.i"
%include "dolfin/mesh/CellType.h"
%include "dolfin/mesh/MeshTopology.h"
%include "dolfin/mesh/MeshGeometry.h"
%include "dolfin/mesh/MeshDomains.h"
%include "dolfin/mesh/MeshData.h"
%include "dolfin/mesh/Mesh.h"
%include "dolfin/mesh/MeshEntity.h"
%include "dolfin/mesh/MeshEntityIterator.h"
%include "dolfin/mesh/MeshEntityIteratorBase.h"
%include "dolfin/mesh/SubsetIterator.h"
%include "dolfin/mesh/Point.h"
%include "dolfin/mesh/Vertex.h"
%include "dolfin/mesh/Edge.h"
%include "dolfin/mesh/Face.h"
%include "dolfin/mesh/Facet.h"
%include "dolfin/mesh/Cell.h"
%include "dolfin/mesh/FacetCell.h"
%include "dolfin/mesh/MeshConnectivity.h"
%include "dolfin/mesh/MeshEditor.h"
%include "dolfin/mesh/DynamicMeshEditor.h"
%include "dolfin/mesh/LocalMeshValueCollection.h"
%include "dolfin/mesh/MeshFunction.h"
%include "dolfin/mesh/MeshPartitioning.h"
%include "dolfin/mesh/MeshValueCollection.h"
%include "dolfin/mesh/MeshColoring.h"
%include "dolfin/mesh/MeshRenumbering.h"
%include "dolfin/mesh/MeshTransformation.h"
%include "dolfin/mesh/LocalMeshData.h"
%include "dolfin/mesh/SubDomain.h"
%include "dolfin/mesh/SubMesh.h"
%include "dolfin/mesh/Restriction.h"
%include "dolfin/mesh/DomainBoundary.h"
%include "dolfin/mesh/BoundaryMesh.h"
%include "dolfin/mesh/PeriodicBoundaryComputation.h"
%include "dolfin/swig/mesh/post.i"

// %include types from submodule generation
%include "dolfin/generation/PolygonalMeshGenerator.h"
%include "dolfin/generation/PolyhedralMeshGenerator.h"
%include "dolfin/generation/Triangulate.h"
%include "dolfin/generation/BoxMesh.h"
%include "dolfin/generation/IntervalMesh.h"
%include "dolfin/generation/Interval.h"
%include "dolfin/generation/RectangleMesh.h"
%include "dolfin/generation/UnitTetrahedronMesh.h"
%include "dolfin/generation/UnitCubeMesh.h"
%include "dolfin/generation/UnitCube.h"
%include "dolfin/generation/UnitIntervalMesh.h"
%include "dolfin/generation/UnitInterval.h"
%include "dolfin/generation/UnitTriangleMesh.h"
%include "dolfin/generation/UnitSquareMesh.h"
%include "dolfin/generation/UnitSquare.h"
%include "dolfin/generation/UnitCircleMesh.h"
%include "dolfin/generation/UnitCircle.h"
%include "dolfin/generation/CSGGeometry.h"
%include "dolfin/generation/CSGMeshGenerator.h"
%include "dolfin/generation/CSGCGALMeshGenerator2D.h"
%include "dolfin/generation/CSGCGALMeshGenerator3D.h"
%include "dolfin/generation/CSGOperators.h"
%include "dolfin/generation/CSGPrimitive.h"
%include "dolfin/generation/CSGPrimitives2D.h"
%include "dolfin/generation/CSGPrimitives3D.h"
%include "dolfin/generation/CSGGeometries3D.h"
%include "dolfin/swig/generation/post.i"

// %include types from submodule refinement
%include "dolfin/refinement/refine.h"

// %include types from submodule graph
%include "dolfin/graph/Graph.h"
%include "dolfin/graph/GraphBuilder.h"
%include "dolfin/graph/BoostGraphOrdering.h"
%include "dolfin/graph/SCOTCH.h"
%include "dolfin/swig/graph/post.i"

// %include types from submodule ale
%include "dolfin/swig/ale/pre.i"
%include "dolfin/ale/ALE.h"
%include "dolfin/ale/MeshDisplacement.h"

