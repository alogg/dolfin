// Automatically generated by FFC, the FEniCS Form Compiler, version 0.1.7.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __CONVECTIONDIFFUSION_H
#define __CONVECTIONDIFFUSION_H

#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace ConvectionDiffusion {

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

    void interpolate(const Function& function, real coefficients[], const AffineMap& map) const
    {
      coefficients[0] = function(map(0.000000000000000e+00, 0.000000000000000e+00));
      coefficients[1] = function(map(1.000000000000000e+00, 0.000000000000000e+00));
      coefficients[2] = function(map(0.000000000000000e+00, 1.000000000000000e+00));
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

    void interpolate(const Function& function, real coefficients[], const AffineMap& map) const
    {
      coefficients[0] = function(map(0.000000000000000e+00, 0.000000000000000e+00));
      coefficients[1] = function(map(1.000000000000000e+00, 0.000000000000000e+00));
      coefficients[2] = function(map(0.000000000000000e+00, 1.000000000000000e+00));
    }

  private:

    unsigned int* tensordims;

  };
    
  class FunctionElement_0 : public dolfin::FiniteElement
  {
  public:

    FunctionElement_0() : dolfin::FiniteElement(), tensordims(0)
    {
      tensordims = new unsigned int [1];
      tensordims[0] = 2;
    }

    ~FunctionElement_0()
    {
      if ( tensordims ) delete [] tensordims;
    }

    inline unsigned int spacedim() const
    {
      return 6;
    }

    inline unsigned int shapedim() const
    {
      return 2;
    }

    inline unsigned int tensordim(unsigned int i) const
    {
      dolfin_assert(i < 1);
      return tensordims[i];
    }

    inline unsigned int rank() const
    {
      return 1;
    }

    void dofmap(int dofs[], const Cell& cell, const Mesh& mesh) const
    {
      dofs[0] = cell.nodeID(0);
      dofs[1] = cell.nodeID(1);
      dofs[2] = cell.nodeID(2);
      int offset = mesh.noNodes();
      dofs[3] = offset + cell.nodeID(0);
      dofs[4] = offset + cell.nodeID(1);
      dofs[5] = offset + cell.nodeID(2);
    }

    void interpolate(const Function& function, real coefficients[], const AffineMap& map) const
    {
      coefficients[0] = function(map(0.000000000000000e+00, 0.000000000000000e+00), 0);
      coefficients[1] = function(map(1.000000000000000e+00, 0.000000000000000e+00), 0);
      coefficients[2] = function(map(0.000000000000000e+00, 1.000000000000000e+00), 0);
      coefficients[3] = function(map(0.000000000000000e+00, 0.000000000000000e+00), 1);
      coefficients[4] = function(map(1.000000000000000e+00, 0.000000000000000e+00), 1);
      coefficients[5] = function(map(0.000000000000000e+00, 1.000000000000000e+00), 1);
    }

  private:

    unsigned int* tensordims;

  };

  BilinearForm(Function& w0, const real& c0, const real& c1) : dolfin::BilinearForm(1), c0(c0), c1(c1)
  {
    // Create finite element for test space
    _test = new TestElement();

    // Create finite element for trial space
    _trial = new TrialElement();
        
    // Add functions
    add(w0, new FunctionElement_0());
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute geometry tensors
    real G0_ = map.det;
    real G1_0_0_0 = map.det*c0*c[0][0]*map.g00;
    real G1_0_0_1 = map.det*c0*c[0][0]*map.g10;
    real G1_0_1_0 = map.det*c0*c[0][1]*map.g00;
    real G1_0_1_1 = map.det*c0*c[0][1]*map.g10;
    real G1_0_2_0 = map.det*c0*c[0][2]*map.g00;
    real G1_0_2_1 = map.det*c0*c[0][2]*map.g10;
    real G1_1_3_0 = map.det*c0*c[0][3]*map.g01;
    real G1_1_3_1 = map.det*c0*c[0][3]*map.g11;
    real G1_1_4_0 = map.det*c0*c[0][4]*map.g01;
    real G1_1_4_1 = map.det*c0*c[0][4]*map.g11;
    real G1_1_5_0 = map.det*c0*c[0][5]*map.g01;
    real G1_1_5_1 = map.det*c0*c[0][5]*map.g11;
    real G2_0_0 = map.det*c0*c1*(map.g00*map.g00 + map.g01*map.g01);
    real G2_0_1 = map.det*c0*c1*(map.g00*map.g10 + map.g01*map.g11);
    real G2_1_0 = map.det*c0*c1*(map.g10*map.g00 + map.g11*map.g01);
    real G2_1_1 = map.det*c0*c1*(map.g10*map.g10 + map.g11*map.g11);

    // Compute element tensor
    block[0] = 8.333333333333318e-02*G0_ - 4.166666666666659e-02*G1_0_0_0 - 4.166666666666659e-02*G1_0_0_1 - 2.083333333333329e-02*G1_0_1_0 - 2.083333333333329e-02*G1_0_1_1 - 2.083333333333329e-02*G1_0_2_0 - 2.083333333333329e-02*G1_0_2_1 - 4.166666666666659e-02*G1_1_3_0 - 4.166666666666659e-02*G1_1_3_1 - 2.083333333333329e-02*G1_1_4_0 - 2.083333333333329e-02*G1_1_4_1 - 2.083333333333329e-02*G1_1_5_0 - 2.083333333333329e-02*G1_1_5_1 + 2.499999999999999e-01*G2_0_0 + 2.499999999999998e-01*G2_0_1 + 2.499999999999998e-01*G2_1_0 + 2.499999999999998e-01*G2_1_1;
    block[1] = 4.166666666666659e-02*G0_ + 4.166666666666659e-02*G1_0_0_0 + 2.083333333333329e-02*G1_0_1_0 + 2.083333333333329e-02*G1_0_2_0 + 4.166666666666659e-02*G1_1_3_0 + 2.083333333333329e-02*G1_1_4_0 + 2.083333333333329e-02*G1_1_5_0 - 2.499999999999999e-01*G2_0_0 - 2.499999999999998e-01*G2_1_0;
    block[2] = 4.166666666666658e-02*G0_ + 4.166666666666659e-02*G1_0_0_1 + 2.083333333333329e-02*G1_0_1_1 + 2.083333333333329e-02*G1_0_2_1 + 4.166666666666659e-02*G1_1_3_1 + 2.083333333333329e-02*G1_1_4_1 + 2.083333333333329e-02*G1_1_5_1 - 2.499999999999998e-01*G2_0_1 - 2.499999999999998e-01*G2_1_1;
    block[3] = 4.166666666666659e-02*G0_ - 2.083333333333330e-02*G1_0_0_0 - 2.083333333333329e-02*G1_0_0_1 - 4.166666666666659e-02*G1_0_1_0 - 4.166666666666659e-02*G1_0_1_1 - 2.083333333333330e-02*G1_0_2_0 - 2.083333333333329e-02*G1_0_2_1 - 2.083333333333330e-02*G1_1_3_0 - 2.083333333333329e-02*G1_1_3_1 - 4.166666666666659e-02*G1_1_4_0 - 4.166666666666659e-02*G1_1_4_1 - 2.083333333333330e-02*G1_1_5_0 - 2.083333333333329e-02*G1_1_5_1 - 2.499999999999999e-01*G2_0_0 - 2.499999999999998e-01*G2_0_1;
    block[4] = 8.333333333333318e-02*G0_ + 2.083333333333330e-02*G1_0_0_0 + 4.166666666666659e-02*G1_0_1_0 + 2.083333333333330e-02*G1_0_2_0 + 2.083333333333330e-02*G1_1_3_0 + 4.166666666666659e-02*G1_1_4_0 + 2.083333333333330e-02*G1_1_5_0 + 2.499999999999999e-01*G2_0_0;
    block[5] = 4.166666666666659e-02*G0_ + 2.083333333333329e-02*G1_0_0_1 + 4.166666666666659e-02*G1_0_1_1 + 2.083333333333329e-02*G1_0_2_1 + 2.083333333333329e-02*G1_1_3_1 + 4.166666666666659e-02*G1_1_4_1 + 2.083333333333329e-02*G1_1_5_1 + 2.499999999999998e-01*G2_0_1;
    block[6] = 4.166666666666658e-02*G0_ - 2.083333333333329e-02*G1_0_0_0 - 2.083333333333329e-02*G1_0_0_1 - 2.083333333333330e-02*G1_0_1_0 - 2.083333333333329e-02*G1_0_1_1 - 4.166666666666659e-02*G1_0_2_0 - 4.166666666666658e-02*G1_0_2_1 - 2.083333333333329e-02*G1_1_3_0 - 2.083333333333329e-02*G1_1_3_1 - 2.083333333333330e-02*G1_1_4_0 - 2.083333333333329e-02*G1_1_4_1 - 4.166666666666659e-02*G1_1_5_0 - 4.166666666666658e-02*G1_1_5_1 - 2.499999999999998e-01*G2_1_0 - 2.499999999999998e-01*G2_1_1;
    block[7] = 4.166666666666659e-02*G0_ + 2.083333333333329e-02*G1_0_0_0 + 2.083333333333330e-02*G1_0_1_0 + 4.166666666666659e-02*G1_0_2_0 + 2.083333333333329e-02*G1_1_3_0 + 2.083333333333330e-02*G1_1_4_0 + 4.166666666666659e-02*G1_1_5_0 + 2.499999999999998e-01*G2_1_0;
    block[8] = 8.333333333333316e-02*G0_ + 2.083333333333329e-02*G1_0_0_1 + 2.083333333333329e-02*G1_0_1_1 + 4.166666666666658e-02*G1_0_2_1 + 2.083333333333329e-02*G1_1_3_1 + 2.083333333333329e-02*G1_1_4_1 + 4.166666666666658e-02*G1_1_5_1 + 2.499999999999998e-01*G2_1_1;

  }
        
private:

  const real& c0;  const real& c1;

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

    void interpolate(const Function& function, real coefficients[], const AffineMap& map) const
    {
      coefficients[0] = function(map(0.000000000000000e+00, 0.000000000000000e+00));
      coefficients[1] = function(map(1.000000000000000e+00, 0.000000000000000e+00));
      coefficients[2] = function(map(0.000000000000000e+00, 1.000000000000000e+00));
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

    void interpolate(const Function& function, real coefficients[], const AffineMap& map) const
    {
      coefficients[0] = function(map(0.000000000000000e+00, 0.000000000000000e+00));
      coefficients[1] = function(map(1.000000000000000e+00, 0.000000000000000e+00));
      coefficients[2] = function(map(0.000000000000000e+00, 1.000000000000000e+00));
    }

  private:

    unsigned int* tensordims;

  };
    
  class FunctionElement_1 : public dolfin::FiniteElement
  {
  public:

    FunctionElement_1() : dolfin::FiniteElement(), tensordims(0)
    {
      tensordims = new unsigned int [1];
      tensordims[0] = 2;
    }

    ~FunctionElement_1()
    {
      if ( tensordims ) delete [] tensordims;
    }

    inline unsigned int spacedim() const
    {
      return 6;
    }

    inline unsigned int shapedim() const
    {
      return 2;
    }

    inline unsigned int tensordim(unsigned int i) const
    {
      dolfin_assert(i < 1);
      return tensordims[i];
    }

    inline unsigned int rank() const
    {
      return 1;
    }

    void dofmap(int dofs[], const Cell& cell, const Mesh& mesh) const
    {
      dofs[0] = cell.nodeID(0);
      dofs[1] = cell.nodeID(1);
      dofs[2] = cell.nodeID(2);
      int offset = mesh.noNodes();
      dofs[3] = offset + cell.nodeID(0);
      dofs[4] = offset + cell.nodeID(1);
      dofs[5] = offset + cell.nodeID(2);
    }

    void interpolate(const Function& function, real coefficients[], const AffineMap& map) const
    {
      coefficients[0] = function(map(0.000000000000000e+00, 0.000000000000000e+00), 0);
      coefficients[1] = function(map(1.000000000000000e+00, 0.000000000000000e+00), 0);
      coefficients[2] = function(map(0.000000000000000e+00, 1.000000000000000e+00), 0);
      coefficients[3] = function(map(0.000000000000000e+00, 0.000000000000000e+00), 1);
      coefficients[4] = function(map(1.000000000000000e+00, 0.000000000000000e+00), 1);
      coefficients[5] = function(map(0.000000000000000e+00, 1.000000000000000e+00), 1);
    }

  private:

    unsigned int* tensordims;

  };
    
  class FunctionElement_2 : public dolfin::FiniteElement
  {
  public:

    FunctionElement_2() : dolfin::FiniteElement(), tensordims(0)
    {
      // Do nothing
    }

    ~FunctionElement_2()
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

    void interpolate(const Function& function, real coefficients[], const AffineMap& map) const
    {
      coefficients[0] = function(map(0.000000000000000e+00, 0.000000000000000e+00));
      coefficients[1] = function(map(1.000000000000000e+00, 0.000000000000000e+00));
      coefficients[2] = function(map(0.000000000000000e+00, 1.000000000000000e+00));
    }

  private:

    unsigned int* tensordims;

  };

  LinearForm(Function& w0, Function& w1, Function& w2, const real& c0, const real& c1) : dolfin::LinearForm(3), c0(c0), c1(c1)
  {
    // Create finite element for test space
    _test = new TestElement();
        
    // Add functions
    add(w0, new FunctionElement_0());
    add(w1, new FunctionElement_1());
    add(w2, new FunctionElement_2());
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute geometry tensors
    real G0_0 = map.det*c[0][0];
    real G0_1 = map.det*c[0][1];
    real G0_2 = map.det*c[0][2];
    real G1_0_0_0_0 = map.det*c0*c[1][0]*c[0][0]*map.g00;
    real G1_0_0_0_1 = map.det*c0*c[1][0]*c[0][0]*map.g10;
    real G1_0_0_1_0 = map.det*c0*c[1][0]*c[0][1]*map.g00;
    real G1_0_0_2_1 = map.det*c0*c[1][0]*c[0][2]*map.g10;
    real G1_0_1_0_0 = map.det*c0*c[1][1]*c[0][0]*map.g00;
    real G1_0_1_0_1 = map.det*c0*c[1][1]*c[0][0]*map.g10;
    real G1_0_1_1_0 = map.det*c0*c[1][1]*c[0][1]*map.g00;
    real G1_0_1_2_1 = map.det*c0*c[1][1]*c[0][2]*map.g10;
    real G1_0_2_0_0 = map.det*c0*c[1][2]*c[0][0]*map.g00;
    real G1_0_2_0_1 = map.det*c0*c[1][2]*c[0][0]*map.g10;
    real G1_0_2_1_0 = map.det*c0*c[1][2]*c[0][1]*map.g00;
    real G1_0_2_2_1 = map.det*c0*c[1][2]*c[0][2]*map.g10;
    real G1_1_3_0_0 = map.det*c0*c[1][3]*c[0][0]*map.g01;
    real G1_1_3_0_1 = map.det*c0*c[1][3]*c[0][0]*map.g11;
    real G1_1_3_1_0 = map.det*c0*c[1][3]*c[0][1]*map.g01;
    real G1_1_3_2_1 = map.det*c0*c[1][3]*c[0][2]*map.g11;
    real G1_1_4_0_0 = map.det*c0*c[1][4]*c[0][0]*map.g01;
    real G1_1_4_0_1 = map.det*c0*c[1][4]*c[0][0]*map.g11;
    real G1_1_4_1_0 = map.det*c0*c[1][4]*c[0][1]*map.g01;
    real G1_1_4_2_1 = map.det*c0*c[1][4]*c[0][2]*map.g11;
    real G1_1_5_0_0 = map.det*c0*c[1][5]*c[0][0]*map.g01;
    real G1_1_5_0_1 = map.det*c0*c[1][5]*c[0][0]*map.g11;
    real G1_1_5_1_0 = map.det*c0*c[1][5]*c[0][1]*map.g01;
    real G1_1_5_2_1 = map.det*c0*c[1][5]*c[0][2]*map.g11;
    real G2_0_0_0 = map.det*c0*c1*c[0][0]*(map.g00*map.g00 + map.g01*map.g01);
    real G2_0_0_1 = map.det*c0*c1*c[0][0]*(map.g00*map.g10 + map.g01*map.g11);
    real G2_0_1_0 = map.det*c0*c1*c[0][1]*(map.g00*map.g00 + map.g01*map.g01);
    real G2_0_2_1 = map.det*c0*c1*c[0][2]*(map.g00*map.g10 + map.g01*map.g11);
    real G2_1_0_0 = map.det*c0*c1*c[0][0]*(map.g10*map.g00 + map.g11*map.g01);
    real G2_1_0_1 = map.det*c0*c1*c[0][0]*(map.g10*map.g10 + map.g11*map.g11);
    real G2_1_1_0 = map.det*c0*c1*c[0][1]*(map.g10*map.g00 + map.g11*map.g01);
    real G2_1_2_1 = map.det*c0*c1*c[0][2]*(map.g10*map.g10 + map.g11*map.g11);
    real G3_0 = map.det*c[2][0];
    real G3_1 = map.det*c[2][1];
    real G3_2 = map.det*c[2][2];

    // Compute element tensor
    block[0] = 8.333333333333318e-02*G0_0 + 4.166666666666659e-02*G0_1 + 4.166666666666658e-02*G0_2 + 4.166666666666659e-02*G1_0_0_0_0 + 4.166666666666659e-02*G1_0_0_0_1 - 4.166666666666659e-02*G1_0_0_1_0 - 4.166666666666659e-02*G1_0_0_2_1 + 2.083333333333329e-02*G1_0_1_0_0 + 2.083333333333329e-02*G1_0_1_0_1 - 2.083333333333329e-02*G1_0_1_1_0 - 2.083333333333329e-02*G1_0_1_2_1 + 2.083333333333329e-02*G1_0_2_0_0 + 2.083333333333329e-02*G1_0_2_0_1 - 2.083333333333329e-02*G1_0_2_1_0 - 2.083333333333329e-02*G1_0_2_2_1 + 4.166666666666659e-02*G1_1_3_0_0 + 4.166666666666659e-02*G1_1_3_0_1 - 4.166666666666659e-02*G1_1_3_1_0 - 4.166666666666659e-02*G1_1_3_2_1 + 2.083333333333329e-02*G1_1_4_0_0 + 2.083333333333329e-02*G1_1_4_0_1 - 2.083333333333329e-02*G1_1_4_1_0 - 2.083333333333329e-02*G1_1_4_2_1 + 2.083333333333329e-02*G1_1_5_0_0 + 2.083333333333329e-02*G1_1_5_0_1 - 2.083333333333329e-02*G1_1_5_1_0 - 2.083333333333329e-02*G1_1_5_2_1 - 2.499999999999999e-01*G2_0_0_0 - 2.499999999999998e-01*G2_0_0_1 + 2.499999999999999e-01*G2_0_1_0 + 2.499999999999998e-01*G2_0_2_1 - 2.499999999999998e-01*G2_1_0_0 - 2.499999999999998e-01*G2_1_0_1 + 2.499999999999998e-01*G2_1_1_0 + 2.499999999999998e-01*G2_1_2_1 + 8.333333333333318e-02*G3_0 + 4.166666666666659e-02*G3_1 + 4.166666666666658e-02*G3_2;
    block[1] = 4.166666666666659e-02*G0_0 + 8.333333333333318e-02*G0_1 + 4.166666666666659e-02*G0_2 + 2.083333333333330e-02*G1_0_0_0_0 + 2.083333333333329e-02*G1_0_0_0_1 - 2.083333333333330e-02*G1_0_0_1_0 - 2.083333333333329e-02*G1_0_0_2_1 + 4.166666666666659e-02*G1_0_1_0_0 + 4.166666666666659e-02*G1_0_1_0_1 - 4.166666666666659e-02*G1_0_1_1_0 - 4.166666666666659e-02*G1_0_1_2_1 + 2.083333333333330e-02*G1_0_2_0_0 + 2.083333333333329e-02*G1_0_2_0_1 - 2.083333333333330e-02*G1_0_2_1_0 - 2.083333333333329e-02*G1_0_2_2_1 + 2.083333333333330e-02*G1_1_3_0_0 + 2.083333333333329e-02*G1_1_3_0_1 - 2.083333333333330e-02*G1_1_3_1_0 - 2.083333333333329e-02*G1_1_3_2_1 + 4.166666666666659e-02*G1_1_4_0_0 + 4.166666666666659e-02*G1_1_4_0_1 - 4.166666666666659e-02*G1_1_4_1_0 - 4.166666666666659e-02*G1_1_4_2_1 + 2.083333333333330e-02*G1_1_5_0_0 + 2.083333333333329e-02*G1_1_5_0_1 - 2.083333333333330e-02*G1_1_5_1_0 - 2.083333333333329e-02*G1_1_5_2_1 + 2.499999999999999e-01*G2_0_0_0 + 2.499999999999998e-01*G2_0_0_1 - 2.499999999999999e-01*G2_0_1_0 - 2.499999999999998e-01*G2_0_2_1 + 4.166666666666659e-02*G3_0 + 8.333333333333318e-02*G3_1 + 4.166666666666659e-02*G3_2;
    block[2] = 4.166666666666658e-02*G0_0 + 4.166666666666659e-02*G0_1 + 8.333333333333316e-02*G0_2 + 2.083333333333329e-02*G1_0_0_0_0 + 2.083333333333329e-02*G1_0_0_0_1 - 2.083333333333329e-02*G1_0_0_1_0 - 2.083333333333329e-02*G1_0_0_2_1 + 2.083333333333330e-02*G1_0_1_0_0 + 2.083333333333329e-02*G1_0_1_0_1 - 2.083333333333330e-02*G1_0_1_1_0 - 2.083333333333329e-02*G1_0_1_2_1 + 4.166666666666659e-02*G1_0_2_0_0 + 4.166666666666658e-02*G1_0_2_0_1 - 4.166666666666659e-02*G1_0_2_1_0 - 4.166666666666658e-02*G1_0_2_2_1 + 2.083333333333329e-02*G1_1_3_0_0 + 2.083333333333329e-02*G1_1_3_0_1 - 2.083333333333329e-02*G1_1_3_1_0 - 2.083333333333329e-02*G1_1_3_2_1 + 2.083333333333330e-02*G1_1_4_0_0 + 2.083333333333329e-02*G1_1_4_0_1 - 2.083333333333330e-02*G1_1_4_1_0 - 2.083333333333329e-02*G1_1_4_2_1 + 4.166666666666659e-02*G1_1_5_0_0 + 4.166666666666658e-02*G1_1_5_0_1 - 4.166666666666659e-02*G1_1_5_1_0 - 4.166666666666658e-02*G1_1_5_2_1 + 2.499999999999998e-01*G2_1_0_0 + 2.499999999999998e-01*G2_1_0_1 - 2.499999999999998e-01*G2_1_1_0 - 2.499999999999998e-01*G2_1_2_1 + 4.166666666666658e-02*G3_0 + 4.166666666666659e-02*G3_1 + 8.333333333333316e-02*G3_2;

  }
        
private:

  const real& c0;  const real& c1;

};

} }

#endif
