// Automatically generated by FFC, the FEniCS Form Compiler, version 0.1.8.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __ELASTICITYUPDATED_H
#define __ELASTICITYUPDATED_H

#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace ElasticityUpdated {

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
    
  class FunctionElement_1 : public dolfin::FiniteElement
  {
  public:

    FunctionElement_1() : dolfin::FiniteElement(), tensordims(0)
    {
      tensordims = new unsigned int [1];
      tensordims[0] = 3;
    }

    ~FunctionElement_1()
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
    
  class FunctionElement_2 : public dolfin::FiniteElement
  {
  public:

    FunctionElement_2() : dolfin::FiniteElement(), tensordims(0)
    {
      tensordims = new unsigned int [1];
      tensordims[0] = 3;
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
    
  class FunctionElement_3 : public dolfin::FiniteElement
  {
  public:

    FunctionElement_3() : dolfin::FiniteElement(), tensordims(0)
    {
      tensordims = new unsigned int [1];
      tensordims[0] = 3;
    }

    ~FunctionElement_3()
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

  LinearForm(Function& w0, Function& w1, Function& w2, Function& w3) : dolfin::LinearForm(4)
  {
    // Create finite element for test space
    _test = new TestElement();
        
    // Add functions
    add(w0, new FunctionElement_0());
    add(w1, new FunctionElement_1());
    add(w2, new FunctionElement_2());
    add(w3, new FunctionElement_3());
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute geometry tensors
    real G0_0 = map.det*c[0][0];
    real G0_1 = map.det*c[0][1];
    real G0_2 = map.det*c[0][2];
    real G0_3 = map.det*c[0][3];
    real G0_4 = map.det*c[0][4];
    real G0_5 = map.det*c[0][5];
    real G0_6 = map.det*c[0][6];
    real G0_7 = map.det*c[0][7];
    real G0_8 = map.det*c[0][8];
    real G0_9 = map.det*c[0][9];
    real G0_10 = map.det*c[0][10];
    real G0_11 = map.det*c[0][11];
    real G1_0_0 = map.det*c[1][0]*map.g00;
    real G1_0_1 = map.det*c[1][0]*map.g10;
    real G1_0_2 = map.det*c[1][0]*map.g20;
    real G2_1_0 = map.det*c[1][1]*map.g01;
    real G2_1_1 = map.det*c[1][1]*map.g11;
    real G2_1_2 = map.det*c[1][1]*map.g21;
    real G3_2_0 = map.det*c[1][2]*map.g02;
    real G3_2_1 = map.det*c[1][2]*map.g12;
    real G3_2_2 = map.det*c[1][2]*map.g22;
    real G4_0_0 = map.det*c[2][0]*map.g00;
    real G4_0_1 = map.det*c[2][0]*map.g10;
    real G4_0_2 = map.det*c[2][0]*map.g20;
    real G5_1_0 = map.det*c[2][1]*map.g01;
    real G5_1_1 = map.det*c[2][1]*map.g11;
    real G5_1_2 = map.det*c[2][1]*map.g21;
    real G6_2_0 = map.det*c[2][2]*map.g02;
    real G6_2_1 = map.det*c[2][2]*map.g12;
    real G6_2_2 = map.det*c[2][2]*map.g22;
    real G7_0_0 = map.det*c[3][0]*map.g00;
    real G7_0_1 = map.det*c[3][0]*map.g10;
    real G7_0_2 = map.det*c[3][0]*map.g20;
    real G8_1_0 = map.det*c[3][1]*map.g01;
    real G8_1_1 = map.det*c[3][1]*map.g11;
    real G8_1_2 = map.det*c[3][1]*map.g21;
    real G9_2_0 = map.det*c[3][2]*map.g02;
    real G9_2_1 = map.det*c[3][2]*map.g12;
    real G9_2_2 = map.det*c[3][2]*map.g22;

    // Compute element tensor
    block[0] = 1.666666666666662e-02*G0_0 + 8.333333333333309e-03*G0_1 + 8.333333333333309e-03*G0_2 + 8.333333333333312e-03*G0_3 + 1.666666666666665e-01*G1_0_0 + 1.666666666666665e-01*G1_0_1 + 1.666666666666664e-01*G1_0_2 + 1.666666666666665e-01*G2_1_0 + 1.666666666666665e-01*G2_1_1 + 1.666666666666664e-01*G2_1_2 + 1.666666666666665e-01*G3_2_0 + 1.666666666666665e-01*G3_2_1 + 1.666666666666664e-01*G3_2_2;
    block[1] = 8.333333333333309e-03*G0_0 + 1.666666666666661e-02*G0_1 + 8.333333333333309e-03*G0_2 + 8.333333333333311e-03*G0_3 - 1.666666666666665e-01*G1_0_0 - 1.666666666666665e-01*G2_1_0 - 1.666666666666665e-01*G3_2_0;
    block[2] = 8.333333333333309e-03*G0_0 + 8.333333333333309e-03*G0_1 + 1.666666666666662e-02*G0_2 + 8.333333333333311e-03*G0_3 - 1.666666666666665e-01*G1_0_1 - 1.666666666666665e-01*G2_1_1 - 1.666666666666665e-01*G3_2_1;
    block[3] = 8.333333333333312e-03*G0_0 + 8.333333333333311e-03*G0_1 + 8.333333333333312e-03*G0_2 + 1.666666666666662e-02*G0_3 - 1.666666666666665e-01*G1_0_2 - 1.666666666666665e-01*G2_1_2 - 1.666666666666665e-01*G3_2_2;
    block[4] = 1.666666666666662e-02*G0_4 + 8.333333333333309e-03*G0_5 + 8.333333333333307e-03*G0_6 + 8.333333333333312e-03*G0_7 + 1.666666666666665e-01*G4_0_0 + 1.666666666666665e-01*G4_0_1 + 1.666666666666665e-01*G4_0_2 + 1.666666666666665e-01*G5_1_0 + 1.666666666666665e-01*G5_1_1 + 1.666666666666665e-01*G5_1_2 + 1.666666666666665e-01*G6_2_0 + 1.666666666666665e-01*G6_2_1 + 1.666666666666665e-01*G6_2_2;
    block[5] = 8.333333333333309e-03*G0_4 + 1.666666666666662e-02*G0_5 + 8.333333333333311e-03*G0_6 + 8.333333333333312e-03*G0_7 - 1.666666666666665e-01*G4_0_0 - 1.666666666666665e-01*G5_1_0 - 1.666666666666665e-01*G6_2_0;
    block[6] = 8.333333333333311e-03*G0_4 + 8.333333333333312e-03*G0_5 + 1.666666666666662e-02*G0_6 + 8.333333333333314e-03*G0_7 - 1.666666666666665e-01*G4_0_1 - 1.666666666666665e-01*G5_1_1 - 1.666666666666665e-01*G6_2_1;
    block[7] = 8.333333333333312e-03*G0_4 + 8.333333333333312e-03*G0_5 + 8.333333333333314e-03*G0_6 + 1.666666666666662e-02*G0_7 - 1.666666666666665e-01*G4_0_2 - 1.666666666666665e-01*G5_1_2 - 1.666666666666665e-01*G6_2_2;
    block[8] = 1.666666666666662e-02*G0_8 + 8.333333333333309e-03*G0_9 + 8.333333333333307e-03*G0_10 + 8.333333333333312e-03*G0_11 + 1.666666666666665e-01*G7_0_0 + 1.666666666666665e-01*G7_0_1 + 1.666666666666665e-01*G7_0_2 + 1.666666666666665e-01*G8_1_0 + 1.666666666666665e-01*G8_1_1 + 1.666666666666665e-01*G8_1_2 + 1.666666666666665e-01*G9_2_0 + 1.666666666666665e-01*G9_2_1 + 1.666666666666665e-01*G9_2_2;
    block[9] = 8.333333333333309e-03*G0_8 + 1.666666666666662e-02*G0_9 + 8.333333333333311e-03*G0_10 + 8.333333333333312e-03*G0_11 - 1.666666666666665e-01*G7_0_0 - 1.666666666666665e-01*G8_1_0 - 1.666666666666665e-01*G9_2_0;
    block[10] = 8.333333333333311e-03*G0_8 + 8.333333333333312e-03*G0_9 + 1.666666666666662e-02*G0_10 + 8.333333333333314e-03*G0_11 - 1.666666666666665e-01*G7_0_1 - 1.666666666666665e-01*G8_1_1 - 1.666666666666665e-01*G9_2_1;
    block[11] = 8.333333333333312e-03*G0_8 + 8.333333333333312e-03*G0_9 + 8.333333333333314e-03*G0_10 + 1.666666666666662e-02*G0_11 - 1.666666666666665e-01*G7_0_2 - 1.666666666666665e-01*G8_1_2 - 1.666666666666665e-01*G9_2_2;

  }

};

} }

#endif
