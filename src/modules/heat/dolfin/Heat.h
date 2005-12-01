// Automatically generated by FFC, the FEniCS Form Compiler, version 0.2.3.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __HEAT_H
#define __HEAT_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/Vector.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace Heat {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
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
    
  class TrialElement : public dolfin::FiniteElement
  {
  public:
  
    TrialElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
    {
      // Element is scalar, don't need to initialize tensordims
  
      // Element is simple, don't need to initialize subelements
    }
  
    ~TrialElement()
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
  
  BilinearForm() : dolfin::BilinearForm(0)
  {
    // Create finite element for test space
    _test = new TestElement();

    // Create finite element for trial space
    _trial = new TrialElement();
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute geometry tensors
    const real G0_ = map.det;

    // Compute element tensor
    block[0] = 8.333333333333318e-02*G0_;
    block[1] = 4.166666666666659e-02*G0_;
    block[2] = 4.166666666666658e-02*G0_;
    block[3] = 4.166666666666659e-02*G0_;
    block[4] = 8.333333333333318e-02*G0_;
    block[5] = 4.166666666666659e-02*G0_;
    block[6] = 4.166666666666658e-02*G0_;
    block[7] = 4.166666666666659e-02*G0_;
    block[8] = 8.333333333333316e-02*G0_;
  }

};

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
    
  class FunctionElement_0 : public dolfin::FiniteElement
  {
  public:
  
    FunctionElement_0() : dolfin::FiniteElement(), tensordims(0), subelements(0)
    {
      // Element is scalar, don't need to initialize tensordims
  
      // Element is simple, don't need to initialize subelements
    }
  
    ~FunctionElement_0()
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
    
  class FunctionElement_1 : public dolfin::FiniteElement
  {
  public:
  
    FunctionElement_1() : dolfin::FiniteElement(), tensordims(0), subelements(0)
    {
      // Element is scalar, don't need to initialize tensordims
  
      // Element is simple, don't need to initialize subelements
    }
  
    ~FunctionElement_1()
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
  
  LinearForm(Function& w0, Function& w1) : dolfin::LinearForm(2)
  {
    // Create finite element for test space
    _test = new TestElement();

    // Add functions
    add(w0, new FunctionElement_0());
    add(w1, new FunctionElement_1());
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute coefficients
    const real c0_0 = c[0][0];
    const real c0_1 = c[0][1];
    const real c0_2 = c[0][2];
    const real c1_0 = c[1][0];
    const real c1_1 = c[1][1];
    const real c1_2 = c[1][2];

    // Compute geometry tensors
    const real G0_0_0_0 = map.det*c0_0*map.g00*map.g00 + map.det*c0_0*map.g01*map.g01;
    const real G0_0_0_1 = map.det*c0_0*map.g00*map.g10 + map.det*c0_0*map.g01*map.g11;
    const real G0_0_1_0 = map.det*c0_0*map.g10*map.g00 + map.det*c0_0*map.g11*map.g01;
    const real G0_0_1_1 = map.det*c0_0*map.g10*map.g10 + map.det*c0_0*map.g11*map.g11;
    const real G0_1_0_0 = map.det*c0_1*map.g00*map.g00 + map.det*c0_1*map.g01*map.g01;
    const real G0_1_0_1 = map.det*c0_1*map.g00*map.g10 + map.det*c0_1*map.g01*map.g11;
    const real G0_2_1_0 = map.det*c0_2*map.g10*map.g00 + map.det*c0_2*map.g11*map.g01;
    const real G0_2_1_1 = map.det*c0_2*map.g10*map.g10 + map.det*c0_2*map.g11*map.g11;
    const real G1_0 = map.det*c1_0;
    const real G1_1 = map.det*c1_1;
    const real G1_2 = map.det*c1_2;

    // Compute element tensor
    block[0] = -4.999999999999998e-01*G0_0_0_0 - 4.999999999999997e-01*G0_0_0_1 - 4.999999999999997e-01*G0_0_1_0 - 4.999999999999996e-01*G0_0_1_1 + 4.999999999999998e-01*G0_1_0_0 + 4.999999999999997e-01*G0_1_0_1 + 4.999999999999997e-01*G0_2_1_0 + 4.999999999999996e-01*G0_2_1_1 + 8.333333333333318e-02*G1_0 + 4.166666666666659e-02*G1_1 + 4.166666666666658e-02*G1_2;
    block[1] = 4.999999999999998e-01*G0_0_0_0 + 4.999999999999997e-01*G0_0_1_0 - 4.999999999999998e-01*G0_1_0_0 - 4.999999999999997e-01*G0_2_1_0 + 4.166666666666659e-02*G1_0 + 8.333333333333318e-02*G1_1 + 4.166666666666659e-02*G1_2;
    block[2] = 4.999999999999997e-01*G0_0_0_1 + 4.999999999999996e-01*G0_0_1_1 - 4.999999999999997e-01*G0_1_0_1 - 4.999999999999996e-01*G0_2_1_1 + 4.166666666666658e-02*G1_0 + 4.166666666666659e-02*G1_1 + 8.333333333333316e-02*G1_2;
  }

};

} }

#endif
