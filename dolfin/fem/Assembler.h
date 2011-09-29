// Copyright (C) 2007-2009 Anders Logg
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
// Modified by Garth N. Wells, 2007-2008.
// Modified by Ola Skavhaug, 2008.
//
// First added:  2007-01-17
// Last changed: 2011-03-11

#ifndef __ASSEMBLER_H
#define __ASSEMBLER_H

#include <vector>
#include <dolfin/common/types.h>

namespace dolfin
{

  // Forward declarations
  class GenericTensor;
  class Form;
  class SubDomain;
  class UFC;
  template<typename T> class MeshFunction;

  /// This class provides automated assembly of linear systems, or
  /// more generally, assembly of a sparse tensor from a given
  /// variational form.
  ///
  /// Subdomains for cells and facets may be specified in a number of
  /// different ways:
  ///
  /// 1. By explicitly passing _MeshFunction_ (as pointers) to the
  ///    assemble functions
  ///
  /// 2. By assigning subdomain indicators specified by _MeshFunction_
  ///    to the _Form_ being assembled:
  ///
  ///    .. code-block:: c++
  ///
  ///        form.dx = cell_domains
  ///        form.ds = exterior_facet_domains
  ///        form.dS = interior_facet_domains
  ///
  /// 3. By _MeshFunction_ stored in _MeshData_ as
  ///
  ///    * "cell_domains"
  ///    * "exterior_facet_domains"
  ///    * "interior_facet_domains"
  ///
  /// 4. By specifying a _SubDomain_ which specifies the domain numbered
  ///    as 0 (with the rest treated as domain number 1)
  ///
  /// Note that (1) overrides (2), which overrides (3).

  class Assembler
  {
  public:

    /// Assemble tensor from given form
    ///
    /// *Arguments*
    ///     A (_GenericTensor_)
    ///         The tensor to assemble.
    ///     a (_Form_)
    ///         The form to assemble the tensor from.
    ///     reset_sparsity (bool)
    ///         Optional argument: Default value is true.
    ///     add_values (bool)
    ///         Optional argument: Default value is false.
    static void assemble(GenericTensor& A,
                         const Form& a,
                         bool reset_sparsity=true,
                         bool add_values=false);

    /// Assemble tensor from given form on subdomain
    ///
    /// *Arguments*
    ///     A (_GenericTensor_)
    ///         The tensor to assemble.
    ///     a (_Form_)
    ///         The form to assemble the tensor from.
    ///     sub_domain (_SubDomain_)
    ///         The subdomain to assemble on.
    ///     reset_sparsity (bool)
    ///         Optional argument: Default value is true.
    ///     add_values (bool)
    ///         Optional argument: Default value is false.
    static void assemble(GenericTensor& A,
                         const Form& a,
                         const SubDomain& sub_domain,
                         bool reset_sparsity=true,
                         bool add_values=false);

    /// Assemble tensor from given form on subdomains
    ///
    /// *Arguments*
    ///     A (_GenericTensor_)
    ///         The tensor to assemble.
    ///     a (_Form_)
    ///         The form to assemble the tensor from.
    ///     cell_domains (_MeshFunction_ <uint>)
    ///         Cell domains.
    ///     exterior_facet_domains (_MeshFunction_ <uint>)
    ///         The exterior facet domains.
    ///     interior_facet_domains (_MeshFunction_ <uint>)
    ///         The interior facet domains.
    ///     reset_sparsity (bool)
    ///         Optional argument: Default value is true.
    ///     add_values (bool)
    ///         Optional argument: Default value is false.
    static void assemble(GenericTensor& A,
                         const Form& a,
                         const MeshFunction<uint>* cell_domains,
                         const MeshFunction<uint>* exterior_facet_domains,
                         const MeshFunction<uint>* interior_facet_domains,
                         bool reset_sparsity=true,
                         bool add_values=false);

    /// Assemble tensor from given form over cells. This function is
    /// provided for users who wish to build a customized assembler.
    static void assemble_cells(GenericTensor& A,
                               const Form& a,
                               UFC& ufc,
                               const MeshFunction<uint>* domains,
                               std::vector<double>* values);

    /// Assemble tensor from given form over exterior facets. This
    /// function is provided for users who wish to build a customized
    /// assembler.
    static void assemble_exterior_facets(GenericTensor& A,
                                         const Form& a,
                                         UFC& ufc,
                                         const MeshFunction<uint>* domains,
                                         std::vector<double>* values);

    /// Assemble tensor from given form over interior facets. This
    /// function is provided for users who wish to build a customized
    /// assembler.
    static void assemble_interior_facets(GenericTensor& A,
                                         const Form& a,
                                         UFC& ufc,
                                         const MeshFunction<uint>* domains,
                                         std::vector<double>* values);

  };

}

#endif
