// Automatically generated by FFC, the FEniCS Form Compiler, version 0.1.8.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __ELASTICITYUPDATEDSIGMA0_H
#define __ELASTICITYUPDATEDSIGMA0_H

#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace ElasticityUpdatedSigma0 {

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
      tensordims = new unsigned int [1];
      tensordims[0] = 3;
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
      return 3;
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
      dofs[0] = 3*cell.id() + 0;
      dofs[1] = 3*cell.id() + 1;
      dofs[2] = 3*cell.id() + 2;
    }

    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
      points[1] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
      points[2] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
      components[0] = 0;
      components[1] = 1;
      components[2] = 2;
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
      tensordims[0] = 3;
    }

    ~FunctionElement_0()
    {
      if ( tensordims ) delete [] tensordims;
    }

    inline unsigned int spacedim() const
    {
      return 12;
    }

    inline unsigned int shapedim() const
    {
      return 3;
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
      dofs[3] = cell.nodeID(3);
      int offset = mesh.noNodes();
      dofs[4] = offset + cell.nodeID(0);
      dofs[5] = offset + cell.nodeID(1);
      dofs[6] = offset + cell.nodeID(2);
      dofs[7] = offset + cell.nodeID(3);
      offset = offset + mesh.noNodes();
      dofs[8] = offset + cell.nodeID(0);
      dofs[9] = offset + cell.nodeID(1);
      dofs[10] = offset + cell.nodeID(2);
      dofs[11] = offset + cell.nodeID(3);
    }

    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
      points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
      points[4] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[5] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[6] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
      points[7] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
      points[8] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[9] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[10] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
      points[11] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
      components[4] = 1;
      components[5] = 1;
      components[6] = 1;
      components[7] = 1;
      components[8] = 2;
      components[9] = 2;
      components[10] = 2;
      components[11] = 2;
    }

  private:

    unsigned int* tensordims;

  };

  LinearForm(Function& w0, const real& c0, const real& c1) : dolfin::LinearForm(1), c0(c0), c1(c1)
  {
    // Create finite element for test space
    _test = new TestElement();
        
    // Add functions
    add(w0, new FunctionElement_0());
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute geometry tensors
    real G0_0_0_0 = map.det*c0*c[0][0]*map.g00;
    real G0_0_0_1 = map.det*c0*c[0][0]*map.g10;
    real G0_0_0_2 = map.det*c0*c[0][0]*map.g20;
    real G0_0_1_0 = map.det*c0*c[0][1]*map.g00;
    real G0_0_2_1 = map.det*c0*c[0][2]*map.g10;
    real G0_0_3_2 = map.det*c0*c[0][3]*map.g20;
    real G0_1_4_0 = map.det*c0*c[0][4]*map.g01;
    real G0_1_4_1 = map.det*c0*c[0][4]*map.g11;
    real G0_1_4_2 = map.det*c0*c[0][4]*map.g21;
    real G0_1_5_0 = map.det*c0*c[0][5]*map.g01;
    real G0_1_6_1 = map.det*c0*c[0][6]*map.g11;
    real G0_1_7_2 = map.det*c0*c[0][7]*map.g21;
    real G0_2_8_0 = map.det*c0*c[0][8]*map.g02;
    real G0_2_8_1 = map.det*c0*c[0][8]*map.g12;
    real G0_2_8_2 = map.det*c0*c[0][8]*map.g22;
    real G0_2_9_0 = map.det*c0*c[0][9]*map.g02;
    real G0_2_10_1 = map.det*c0*c[0][10]*map.g12;
    real G0_2_11_2 = map.det*c0*c[0][11]*map.g22;
    real G1_0_0 = map.det*c1*c[0][0]*map.g00;
    real G1_0_1 = map.det*c1*c[0][0]*map.g10;
    real G1_0_2 = map.det*c1*c[0][0]*map.g20;
    real G1_1_0 = map.det*c1*c[0][1]*map.g00;
    real G1_2_1 = map.det*c1*c[0][2]*map.g10;
    real G1_3_2 = map.det*c1*c[0][3]*map.g20;
    real G1_4_0 = map.det*c1*c[0][4]*map.g00;
    real G1_4_1 = map.det*c1*c[0][4]*map.g10;
    real G1_4_2 = map.det*c1*c[0][4]*map.g20;
    real G1_5_0 = map.det*c1*c[0][5]*map.g00;
    real G1_6_1 = map.det*c1*c[0][6]*map.g10;
    real G1_7_2 = map.det*c1*c[0][7]*map.g20;
    real G1_8_0 = map.det*c1*c[0][8]*map.g00;
    real G1_8_1 = map.det*c1*c[0][8]*map.g10;
    real G1_8_2 = map.det*c1*c[0][8]*map.g20;
    real G1_9_0 = map.det*c1*c[0][9]*map.g00;
    real G1_10_1 = map.det*c1*c[0][10]*map.g10;
    real G1_11_2 = map.det*c1*c[0][11]*map.g20;
    real G2_0_4_0 = map.det*c[0][4]*map.g00;
    real G2_0_4_1 = map.det*c[0][4]*map.g10;
    real G2_0_4_2 = map.det*c[0][4]*map.g20;
    real G2_0_5_0 = map.det*c[0][5]*map.g00;
    real G2_0_6_1 = map.det*c[0][6]*map.g10;
    real G2_0_7_2 = map.det*c[0][7]*map.g20;
    real G2_1_4_0 = map.det*c[0][4]*map.g01;
    real G2_1_4_1 = map.det*c[0][4]*map.g11;
    real G2_1_4_2 = map.det*c[0][4]*map.g21;
    real G2_1_5_0 = map.det*c[0][5]*map.g01;
    real G2_1_6_1 = map.det*c[0][6]*map.g11;
    real G2_1_7_2 = map.det*c[0][7]*map.g21;
    real G2_2_4_0 = map.det*c[0][4]*map.g02;
    real G2_2_4_1 = map.det*c[0][4]*map.g12;
    real G2_2_4_2 = map.det*c[0][4]*map.g22;
    real G2_2_5_0 = map.det*c[0][5]*map.g02;
    real G2_2_6_1 = map.det*c[0][6]*map.g12;
    real G2_2_7_2 = map.det*c[0][7]*map.g22;

    // Compute element tensor
    block[0] = -1.666666666666665e-01*G1_0_0 - 1.666666666666665e-01*G1_0_1 - 1.666666666666664e-01*G1_0_2 + 1.666666666666665e-01*G1_1_0 + 1.666666666666665e-01*G1_2_1 + 1.666666666666665e-01*G1_3_2 - 1.666666666666665e-01*G2_0_4_0 - 1.666666666666665e-01*G2_0_4_1 - 1.666666666666665e-01*G2_0_4_2 + 1.666666666666665e-01*G2_0_5_0 + 1.666666666666665e-01*G2_0_6_1 + 1.666666666666665e-01*G2_0_7_2;
    block[1] = -1.666666666666665e-01*G0_0_0_0 - 1.666666666666665e-01*G0_0_0_1 - 1.666666666666664e-01*G0_0_0_2 + 1.666666666666665e-01*G0_0_1_0 + 1.666666666666665e-01*G0_0_2_1 + 1.666666666666665e-01*G0_0_3_2 - 1.666666666666665e-01*G0_1_4_0 - 1.666666666666665e-01*G0_1_4_1 - 1.666666666666665e-01*G0_1_4_2 + 1.666666666666665e-01*G0_1_5_0 + 1.666666666666665e-01*G0_1_6_1 + 1.666666666666665e-01*G0_1_7_2 - 1.666666666666665e-01*G0_2_8_0 - 1.666666666666665e-01*G0_2_8_1 - 1.666666666666665e-01*G0_2_8_2 + 1.666666666666665e-01*G0_2_9_0 + 1.666666666666665e-01*G0_2_10_1 + 1.666666666666665e-01*G0_2_11_2 - 1.666666666666665e-01*G1_4_0 - 1.666666666666665e-01*G1_4_1 - 1.666666666666665e-01*G1_4_2 + 1.666666666666665e-01*G1_5_0 + 1.666666666666665e-01*G1_6_1 + 1.666666666666665e-01*G1_7_2 - 1.666666666666665e-01*G2_1_4_0 - 1.666666666666665e-01*G2_1_4_1 - 1.666666666666665e-01*G2_1_4_2 + 1.666666666666665e-01*G2_1_5_0 + 1.666666666666665e-01*G2_1_6_1 + 1.666666666666665e-01*G2_1_7_2;
    block[2] = -1.666666666666665e-01*G1_8_0 - 1.666666666666665e-01*G1_8_1 - 1.666666666666665e-01*G1_8_2 + 1.666666666666665e-01*G1_9_0 + 1.666666666666665e-01*G1_10_1 + 1.666666666666665e-01*G1_11_2 - 1.666666666666665e-01*G2_2_4_0 - 1.666666666666665e-01*G2_2_4_1 - 1.666666666666665e-01*G2_2_4_2 + 1.666666666666665e-01*G2_2_5_0 + 1.666666666666665e-01*G2_2_6_1 + 1.666666666666665e-01*G2_2_7_2;

  }
        
private:

  const real& c0;  const real& c1;

};

} }

#endif
