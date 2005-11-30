// Copyright (C) 2005 Anders Logg.
// Licensed under the GNU GPL Version 2.
//
// First added:  2005-05-02
// Last changed: 2005-11-29

#ifndef __FINITE_ELEMENT_H
#define __FINITE_ELEMENT_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Node.h>
#include <dolfin/Point.h>
#include <dolfin/Vector.h>
#include <dolfin/AffineMap.h>

namespace dolfin
{
  /// This is the base class for finite elements automatically
  /// generated by the FEniCS Form Compiler FFC.

  class FiniteElement
  {
  public:
   
    /// Constructor
    FiniteElement();

    /// Destructor
    virtual ~FiniteElement();
    
    /// Return dimension of the finite element space
    virtual unsigned int spacedim() const = 0;

    /// Return dimension of the underlying shape
    virtual unsigned int shapedim() const = 0;

    /// Return vector dimension of the finite element space
    virtual unsigned int tensordim(unsigned int i) const = 0;

    /// Return dimension of (mixed) element (number of sub elements)
    virtual unsigned int elementdim() const = 0;

    /// Return vector dimension of the finite element space
    virtual unsigned int rank() const = 0;
    
    /// Compute map from local to global degrees of freedom
    virtual void dofmap(int dofs[], const Cell& cell, const Mesh& mesh) const = 0;

    /// Compute map from local to global coordinates
    virtual void pointmap(Point points[], unsigned int components[],
			  const AffineMap& map) const = 0;

    // FIXME: Should not have vector here!
    /// Compute map from (vertex, component) to function value
    virtual void vertexeval(real values[], unsigned int vertex, 
			    const Vector& x, const Mesh& mesh) const = 0;

    /// Return given sub element of (mixed) element
    virtual const FiniteElement& operator[] (unsigned int i) const = 0;

    /// Display finite element data
    void disp() const;

  };

}

#endif
