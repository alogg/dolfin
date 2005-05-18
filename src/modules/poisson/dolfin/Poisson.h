// Automatically generated by FFC, the FEniCS Form Compiler, version 0.1.8.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __POISSON_H
#define __POISSON_H

#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace Poisson {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
{
public:
    
  class TestElement : public dolfin::FiniteElement
  {
  public:

    TestElement() : dolfin::FiniteElement(), tensordims(0)
    {
      // Do nothing
    }

    ~TestElement()
    {
      if ( tensordims ) delete [] tensordims;
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

  private:

    unsigned int* tensordims;

  };
    
  class TrialElement : public dolfin::FiniteElement
  {
  public:

    TrialElement() : dolfin::FiniteElement(), tensordims(0)
    {
      // Do nothing
    }

    ~TrialElement()
    {
      if ( tensordims ) delete [] tensordims;
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

  private:

    unsigned int* tensordims;

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
    real G0_0_0 = map.det*(map.g00*map.g00 + map.g01*map.g01);
    real G0_0_1 = map.det*(map.g00*map.g10 + map.g01*map.g11);
    real G0_1_0 = map.det*(map.g10*map.g00 + map.g11*map.g01);
    real G0_1_1 = map.det*(map.g10*map.g10 + map.g11*map.g11);

    // Compute element tensor
    block[0] = 4.999999999999998e-01*G0_0_0 + 4.999999999999997e-01*G0_0_1 + 4.999999999999997e-01*G0_1_0 + 4.999999999999996e-01*G0_1_1;
    block[1] = -4.999999999999998e-01*G0_0_0 - 4.999999999999997e-01*G0_1_0;
    block[2] = -4.999999999999997e-01*G0_0_1 - 4.999999999999996e-01*G0_1_1;
    block[3] = -4.999999999999998e-01*G0_0_0 - 4.999999999999997e-01*G0_0_1;
    block[4] = 4.999999999999998e-01*G0_0_0;
    block[5] = 4.999999999999997e-01*G0_0_1;
    block[6] = -4.999999999999997e-01*G0_1_0 - 4.999999999999996e-01*G0_1_1;
    block[7] = 4.999999999999997e-01*G0_1_0;
    block[8] = 4.999999999999996e-01*G0_1_1;

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

    TestElement() : dolfin::FiniteElement(), tensordims(0)
    {
      // Do nothing
    }

    ~TestElement()
    {
      if ( tensordims ) delete [] tensordims;
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

  private:

    unsigned int* tensordims;

  };
    
  class FunctionElement_0 : public dolfin::FiniteElement
  {
  public:

    FunctionElement_0() : dolfin::FiniteElement(), tensordims(0)
    {
      // Do nothing
    }

    ~FunctionElement_0()
    {
      if ( tensordims ) delete [] tensordims;
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

  private:

    unsigned int* tensordims;

  };

  LinearForm(Function& w0) : dolfin::LinearForm(1)
  {
    // Create finite element for test space
    _test = new TestElement();
        
    // Add functions
    add(w0, new FunctionElement_0());
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute geometry tensors
    real G0_0 = map.det*c[0][0];
    real G0_1 = map.det*c[0][1];
    real G0_2 = map.det*c[0][2];

    // Compute element tensor
    block[0] = 8.333333333333318e-02*G0_0 + 4.166666666666659e-02*G0_1 + 4.166666666666658e-02*G0_2;
    block[1] = 4.166666666666659e-02*G0_0 + 8.333333333333318e-02*G0_1 + 4.166666666666659e-02*G0_2;
    block[2] = 4.166666666666658e-02*G0_0 + 4.166666666666659e-02*G0_1 + 8.333333333333316e-02*G0_2;

  }

};

} }

#endif
