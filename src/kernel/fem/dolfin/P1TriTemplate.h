// Automatically generated by FFC, the FEniCS Form Compiler, version 0.2.3.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __P1TRITEMPLATE_H
#define __P1TRITEMPLATE_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/Vector.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace P1TriTemplate {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class LinearForm : public dolfin::LinearForm
{
public:
  
  class TestElement : public dolfin::FiniteElement
  {
  public:
  
    TestElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
    {
      // Element is scalar, don't need to initialize tensordims
  
      // Element is simple, don't need to initialize subelements
    }
  
    ~TestElement()
    {
      if ( tensordims ) delete [] tensordims;
      if ( subelements )
      {
        for (unsigned int i = 0; i < elementdim(); i++)
          delete subelements[i];
        delete [] subelements;
      }
    }
  
    inline unsigned int spacedim() const
    {
      return 3;
    }
  
    inline unsigned int shapedim() const
    {
      return 2;
    }
  
    inline unsigned int tensordim(unsigned int i) const
    {
      dolfin_error("Element is scalar.");
      return 0;
    }
  
    inline unsigned int elementdim() const
    {
      return 1;
    }
  
    inline unsigned int rank() const
    {
      return 0;
    }
  
    void dofmap(int dofs[], const Cell& cell, const Mesh& mesh) const
    {
      dofs[0] = cell.nodeID(0);
      dofs[1] = cell.nodeID(1);
      dofs[2] = cell.nodeID(2);
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
    }
  
    void vertexeval(real values[], unsigned int vertex, const real x[], const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      values[0] = x[vertex];
    }
  
    const FiniteElement& operator[] (unsigned int i) const
    {
      return *this;
    }
  
    FiniteElement& operator[] (unsigned int i)
    {
      return *this;
    }
    
  private:
  
    unsigned int* tensordims;
    FiniteElement** subelements;
  
  };
  
  LinearForm() : dolfin::LinearForm(0)
  {
    // Create finite element for test space
    _test = new TestElement();
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute geometry tensors
    const real G0_ = map.det;

    // Compute element tensor
    block[0] = 1.666666666666665e-01*G0_;
    block[1] = 1.666666666666665e-01*G0_;
    block[2] = 1.666666666666665e-01*G0_;
  }

};

} }

#endif
