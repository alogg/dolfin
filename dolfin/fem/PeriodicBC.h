// Copyright (C) 2007-2008 Anders Logg.
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Garth N. Wells 2007
// Modified by Johan Hake 2009
//
// First added:  2007-07-08
// Last changed: 2009-09-16

#ifndef __PERIODIC_BC_H
#define __PERIODIC_BC_H

#include <boost/shared_ptr.hpp>
#include <dolfin/common/types.h>
#include "BoundaryCondition.h"

namespace dolfin
{

  class DofMap;
  class Mesh;
  class SubDomain;
  class Form;
  class GenericMatrix;
  class GenericVector;

  /// This class specifies the interface for setting periodic boundary
  /// conditions for partial differential equations,
  ///
  ///    u(x) = u(F^{-1}(x)) on G,
  ///    u(x) = u(F(x))      on H,
  ///
  /// where F : H --> G is a map from a subdomain H to a subdomain G.
  ///
  /// A periodic boundary condition must be defined by the domain G
  /// and the map F pulling coordinates back from H to G. The domain
  /// and the map are both defined by a subclass of SubDomain which
  /// must must overload both the inside() function, which specifies
  /// the points of G, and the map() function, which specifies the map
  /// from the points of H to the points of G.
  ///
  /// The implementation is based on matching degrees of freedom on G
  /// with degrees of freedom on H and only works when the mapping F
  /// is bijective between the sets of coordinates associated with the
  /// two domains. In other words, the nodes (degrees of freedom) must
  /// be aligned on G and H.
  ///
  /// For mixed systems (vector-valued and mixed elements), an
  /// optional set of parameters may be used to specify for which sub
  /// system the boundary condition should be specified.

  class PeriodicBC : public BoundaryCondition
  {
  public:

    /// Create periodic boundary condition for sub domain
    PeriodicBC(const FunctionSpace& V,
               const SubDomain& sub_domain);

    /// Create periodic boundary condition for sub domain
    PeriodicBC(boost::shared_ptr<const FunctionSpace> V,
               boost::shared_ptr<const SubDomain> sub_domain);

    /// Destructor
    ~PeriodicBC();

    /// Apply boundary condition to a matrix
    void apply(GenericMatrix& A) const;

    /// Apply boundary condition to a vector
    void apply(GenericVector& b) const;

    /// Apply boundary condition to a linear system
    void apply(GenericMatrix& A, GenericVector& b) const;

    /// Apply boundary condition to a vector for a nonlinear problem
    void apply(GenericVector& b, const GenericVector& x) const;

    /// Apply boundary condition to a linear system for a nonlinear problem
    void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x) const;

    /// Rebuild mapping between dofs
    void rebuild();

  private:

    // The subdomain
    boost::shared_ptr<const SubDomain> sub_domain;

    // Number of dof pairs
    uint num_dof_pairs;

    // Array of master dofs (size num_dof_pairs)
    uint* master_dofs;

    // Array of slave dofs (size num_dof_pairs)
    uint* slave_dofs;

    // Zeros, used for zeroing entries in right-hand side (size num_dof_pairs)
    double* zeros;

  };

}

#endif
