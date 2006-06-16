// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.2.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __CONVECTIONDIFFUSION_H
#define __CONVECTIONDIFFUSION_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/FiniteElementSpec.h>
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
  
    void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
    {
      nodes[0] = cell.vertexID(0);
      nodes[1] = cell.vertexID(1);
      nodes[2] = cell.vertexID(2);
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
  
    void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      vertex_nodes[0] = vertex;
    }
  
    const FiniteElement& operator[] (unsigned int i) const
    {
      return *this;
    }
  
    FiniteElement& operator[] (unsigned int i)
    {
      return *this;
    }
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 1);
      return s;
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
  
    void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
    {
      nodes[0] = cell.vertexID(0);
      nodes[1] = cell.vertexID(1);
      nodes[2] = cell.vertexID(2);
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
  
    void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      vertex_nodes[0] = vertex;
    }
  
    const FiniteElement& operator[] (unsigned int i) const
    {
      return *this;
    }
  
    FiniteElement& operator[] (unsigned int i)
    {
      return *this;
    }
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 1);
      return s;
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
      tensordims = new unsigned int [1];
      tensordims[0] = 2;
  
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
      return 12;
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
  
    inline unsigned int elementdim() const
    {
      return 1;
    }
  
    inline unsigned int rank() const
    {
      return 1;
    }
  
    void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
    {
      nodes[0] = cell.vertexID(0);
      nodes[1] = cell.vertexID(1);
      nodes[2] = cell.vertexID(2);
      int offset = mesh.numVertices();
      nodes[3] = offset + cell.edgeID(0);
      nodes[4] = offset + cell.edgeID(1);
      nodes[5] = offset + cell.edgeID(2);
      offset = offset + mesh.numEdges();
      nodes[6] = offset + cell.vertexID(0);
      nodes[7] = offset + cell.vertexID(1);
      nodes[8] = offset + cell.vertexID(2);
      offset = offset + mesh.numVertices();
      nodes[9] = offset + cell.edgeID(0);
      nodes[10] = offset + cell.edgeID(1);
      nodes[11] = offset + cell.edgeID(2);
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
      points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
      points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
      points[6] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[7] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[8] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[9] = map(5.000000000000000e-01, 5.000000000000000e-01);
      points[10] = map(0.000000000000000e+00, 5.000000000000000e-01);
      points[11] = map(5.000000000000000e-01, 0.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
      components[4] = 0;
      components[5] = 0;
      components[6] = 1;
      components[7] = 1;
      components[8] = 1;
      components[9] = 1;
      components[10] = 1;
      components[11] = 1;
    }
  
    void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      vertex_nodes[0] = vertex;
      int offset = mesh.numVertices() + mesh.numEdges();
      vertex_nodes[1] = offset + vertex;
    }
  
    const FiniteElement& operator[] (unsigned int i) const
    {
      return *this;
    }
  
    FiniteElement& operator[] (unsigned int i)
    {
      return *this;
    }
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Vector Lagrange", "triangle", 2, 2);
      return s;
    }
    
  private:
  
    unsigned int* tensordims;
    FiniteElement** subelements;
  
  };
  
  BilinearForm(Function& w0) : dolfin::BilinearForm(1)
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
    // Compute coefficients
    const real c0_0 = c[0][0];
    const real c0_1 = c[0][1];
    const real c0_2 = c[0][2];
    const real c0_3 = c[0][3];
    const real c0_4 = c[0][4];
    const real c0_5 = c[0][5];
    const real c0_6 = c[0][6];
    const real c0_7 = c[0][7];
    const real c0_8 = c[0][8];
    const real c0_9 = c[0][9];
    const real c0_10 = c[0][10];
    const real c0_11 = c[0][11];

    // Compute geometry tensors
    const real G0_ = map.det;
    const real G1_0_0 = map.det*c0_0*map.g00;
    const real G1_0_1 = map.det*c0_0*map.g10;
    const real G1_1_0 = map.det*c0_1*map.g00;
    const real G1_1_1 = map.det*c0_1*map.g10;
    const real G1_2_0 = map.det*c0_2*map.g00;
    const real G1_2_1 = map.det*c0_2*map.g10;
    const real G1_3_0 = map.det*c0_3*map.g00;
    const real G1_3_1 = map.det*c0_3*map.g10;
    const real G1_4_0 = map.det*c0_4*map.g00;
    const real G1_4_1 = map.det*c0_4*map.g10;
    const real G1_5_0 = map.det*c0_5*map.g00;
    const real G1_5_1 = map.det*c0_5*map.g10;
    const real G2_6_0 = map.det*c0_6*map.g01;
    const real G2_6_1 = map.det*c0_6*map.g11;
    const real G2_7_0 = map.det*c0_7*map.g01;
    const real G2_7_1 = map.det*c0_7*map.g11;
    const real G2_8_0 = map.det*c0_8*map.g01;
    const real G2_8_1 = map.det*c0_8*map.g11;
    const real G2_9_0 = map.det*c0_9*map.g01;
    const real G2_9_1 = map.det*c0_9*map.g11;
    const real G2_10_0 = map.det*c0_10*map.g01;
    const real G2_10_1 = map.det*c0_10*map.g11;
    const real G2_11_0 = map.det*c0_11*map.g01;
    const real G2_11_1 = map.det*c0_11*map.g11;
    const real G3_0_0 = map.det*map.g00*map.g00 + map.det*map.g01*map.g01;
    const real G3_0_1 = map.det*map.g00*map.g10 + map.det*map.g01*map.g11;
    const real G3_1_0 = map.det*map.g10*map.g00 + map.det*map.g11*map.g01;
    const real G3_1_1 = map.det*map.g10*map.g10 + map.det*map.g11*map.g11;

    // Compute element tensor
    block[0] = 8.333333333333318e-02*G0_ - 4.166666666666660e-04*G1_0_0 - 4.166666666666659e-04*G1_0_1 + 2.083333333333329e-04*G1_1_0 + 2.083333333333329e-04*G1_1_1 + 2.083333333333329e-04*G1_2_0 + 2.083333333333329e-04*G1_2_1 - 8.333333333333318e-04*G1_3_0 - 8.333333333333315e-04*G1_3_1 - 1.666666666666664e-03*G1_4_0 - 1.666666666666664e-03*G1_4_1 - 1.666666666666664e-03*G1_5_0 - 1.666666666666663e-03*G1_5_1 - 4.166666666666662e-04*G2_6_0 - 4.166666666666661e-04*G2_6_1 + 2.083333333333328e-04*G2_7_0 + 2.083333333333328e-04*G2_7_1 + 2.083333333333329e-04*G2_8_0 + 2.083333333333329e-04*G2_8_1 - 8.333333333333318e-04*G2_9_0 - 8.333333333333315e-04*G2_9_1 - 1.666666666666664e-03*G2_10_0 - 1.666666666666664e-03*G2_10_1 - 1.666666666666664e-03*G2_11_0 - 1.666666666666663e-03*G2_11_1 + 6.249999999999997e-05*G3_0_0 + 6.249999999999996e-05*G3_0_1 + 6.249999999999996e-05*G3_1_0 + 6.249999999999995e-05*G3_1_1;
    block[1] = 4.166666666666659e-02*G0_ + 4.166666666666660e-04*G1_0_0 - 2.083333333333329e-04*G1_1_0 - 2.083333333333329e-04*G1_2_0 + 8.333333333333318e-04*G1_3_0 + 1.666666666666664e-03*G1_4_0 + 1.666666666666664e-03*G1_5_0 + 4.166666666666662e-04*G2_6_0 - 2.083333333333328e-04*G2_7_0 - 2.083333333333329e-04*G2_8_0 + 8.333333333333318e-04*G2_9_0 + 1.666666666666664e-03*G2_10_0 + 1.666666666666664e-03*G2_11_0 - 6.249999999999997e-05*G3_0_0 - 6.249999999999996e-05*G3_1_0;
    block[2] = 4.166666666666658e-02*G0_ + 4.166666666666659e-04*G1_0_1 - 2.083333333333329e-04*G1_1_1 - 2.083333333333329e-04*G1_2_1 + 8.333333333333315e-04*G1_3_1 + 1.666666666666664e-03*G1_4_1 + 1.666666666666663e-03*G1_5_1 + 4.166666666666661e-04*G2_6_1 - 2.083333333333328e-04*G2_7_1 - 2.083333333333329e-04*G2_8_1 + 8.333333333333315e-04*G2_9_1 + 1.666666666666664e-03*G2_10_1 + 1.666666666666663e-03*G2_11_1 - 6.249999999999996e-05*G3_0_1 - 6.249999999999995e-05*G3_1_1;
    block[3] = 4.166666666666659e-02*G0_ + 2.083333333333330e-04*G1_0_0 + 2.083333333333329e-04*G1_0_1 - 4.166666666666660e-04*G1_1_0 - 4.166666666666659e-04*G1_1_1 + 2.083333333333328e-04*G1_2_0 + 2.083333333333328e-04*G1_2_1 - 1.666666666666664e-03*G1_3_0 - 1.666666666666664e-03*G1_3_1 - 8.333333333333320e-04*G1_4_0 - 8.333333333333318e-04*G1_4_1 - 1.666666666666664e-03*G1_5_0 - 1.666666666666663e-03*G1_5_1 + 2.083333333333326e-04*G2_6_0 + 2.083333333333325e-04*G2_6_1 - 4.166666666666660e-04*G2_7_0 - 4.166666666666659e-04*G2_7_1 + 2.083333333333328e-04*G2_8_0 + 2.083333333333328e-04*G2_8_1 - 1.666666666666664e-03*G2_9_0 - 1.666666666666664e-03*G2_9_1 - 8.333333333333320e-04*G2_10_0 - 8.333333333333318e-04*G2_10_1 - 1.666666666666664e-03*G2_11_0 - 1.666666666666663e-03*G2_11_1 - 6.249999999999997e-05*G3_0_0 - 6.249999999999996e-05*G3_0_1;
    block[4] = 8.333333333333318e-02*G0_ - 2.083333333333330e-04*G1_0_0 + 4.166666666666660e-04*G1_1_0 - 2.083333333333328e-04*G1_2_0 + 1.666666666666664e-03*G1_3_0 + 8.333333333333320e-04*G1_4_0 + 1.666666666666664e-03*G1_5_0 - 2.083333333333326e-04*G2_6_0 + 4.166666666666660e-04*G2_7_0 - 2.083333333333328e-04*G2_8_0 + 1.666666666666664e-03*G2_9_0 + 8.333333333333320e-04*G2_10_0 + 1.666666666666664e-03*G2_11_0 + 6.249999999999997e-05*G3_0_0;
    block[5] = 4.166666666666659e-02*G0_ - 2.083333333333329e-04*G1_0_1 + 4.166666666666659e-04*G1_1_1 - 2.083333333333328e-04*G1_2_1 + 1.666666666666664e-03*G1_3_1 + 8.333333333333318e-04*G1_4_1 + 1.666666666666663e-03*G1_5_1 - 2.083333333333325e-04*G2_6_1 + 4.166666666666659e-04*G2_7_1 - 2.083333333333328e-04*G2_8_1 + 1.666666666666664e-03*G2_9_1 + 8.333333333333318e-04*G2_10_1 + 1.666666666666663e-03*G2_11_1 + 6.249999999999996e-05*G3_0_1;
    block[6] = 4.166666666666658e-02*G0_ + 2.083333333333330e-04*G1_0_0 + 2.083333333333330e-04*G1_0_1 + 2.083333333333330e-04*G1_1_0 + 2.083333333333329e-04*G1_1_1 - 4.166666666666659e-04*G1_2_0 - 4.166666666666658e-04*G1_2_1 - 1.666666666666664e-03*G1_3_0 - 1.666666666666664e-03*G1_3_1 - 1.666666666666664e-03*G1_4_0 - 1.666666666666664e-03*G1_4_1 - 8.333333333333318e-04*G1_5_0 - 8.333333333333315e-04*G1_5_1 + 2.083333333333327e-04*G2_6_0 + 2.083333333333327e-04*G2_6_1 + 2.083333333333329e-04*G2_7_0 + 2.083333333333329e-04*G2_7_1 - 4.166666666666659e-04*G2_8_0 - 4.166666666666658e-04*G2_8_1 - 1.666666666666664e-03*G2_9_0 - 1.666666666666664e-03*G2_9_1 - 1.666666666666664e-03*G2_10_0 - 1.666666666666664e-03*G2_10_1 - 8.333333333333318e-04*G2_11_0 - 8.333333333333315e-04*G2_11_1 - 6.249999999999996e-05*G3_1_0 - 6.249999999999995e-05*G3_1_1;
    block[7] = 4.166666666666659e-02*G0_ - 2.083333333333330e-04*G1_0_0 - 2.083333333333330e-04*G1_1_0 + 4.166666666666659e-04*G1_2_0 + 1.666666666666664e-03*G1_3_0 + 1.666666666666664e-03*G1_4_0 + 8.333333333333318e-04*G1_5_0 - 2.083333333333327e-04*G2_6_0 - 2.083333333333329e-04*G2_7_0 + 4.166666666666659e-04*G2_8_0 + 1.666666666666664e-03*G2_9_0 + 1.666666666666664e-03*G2_10_0 + 8.333333333333318e-04*G2_11_0 + 6.249999999999996e-05*G3_1_0;
    block[8] = 8.333333333333316e-02*G0_ - 2.083333333333330e-04*G1_0_1 - 2.083333333333329e-04*G1_1_1 + 4.166666666666658e-04*G1_2_1 + 1.666666666666664e-03*G1_3_1 + 1.666666666666664e-03*G1_4_1 + 8.333333333333315e-04*G1_5_1 - 2.083333333333327e-04*G2_6_1 - 2.083333333333329e-04*G2_7_1 + 4.166666666666658e-04*G2_8_1 + 1.666666666666664e-03*G2_9_1 + 1.666666666666664e-03*G2_10_1 + 8.333333333333315e-04*G2_11_1 + 6.249999999999995e-05*G3_1_1;
  }

  // No contribution from the boundary
  void eval(real block[], const AffineMap& map, unsigned int facet) const {}   

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
  
    void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
    {
      nodes[0] = cell.vertexID(0);
      nodes[1] = cell.vertexID(1);
      nodes[2] = cell.vertexID(2);
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
  
    void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      vertex_nodes[0] = vertex;
    }
  
    const FiniteElement& operator[] (unsigned int i) const
    {
      return *this;
    }
  
    FiniteElement& operator[] (unsigned int i)
    {
      return *this;
    }
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 1);
      return s;
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
  
    void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
    {
      nodes[0] = cell.vertexID(0);
      nodes[1] = cell.vertexID(1);
      nodes[2] = cell.vertexID(2);
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
  
    void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      vertex_nodes[0] = vertex;
    }
  
    const FiniteElement& operator[] (unsigned int i) const
    {
      return *this;
    }
  
    FiniteElement& operator[] (unsigned int i)
    {
      return *this;
    }
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 1);
      return s;
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
      tensordims = new unsigned int [1];
      tensordims[0] = 2;
  
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
      return 12;
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
  
    inline unsigned int elementdim() const
    {
      return 1;
    }
  
    inline unsigned int rank() const
    {
      return 1;
    }
  
    void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
    {
      nodes[0] = cell.vertexID(0);
      nodes[1] = cell.vertexID(1);
      nodes[2] = cell.vertexID(2);
      int offset = mesh.numVertices();
      nodes[3] = offset + cell.edgeID(0);
      nodes[4] = offset + cell.edgeID(1);
      nodes[5] = offset + cell.edgeID(2);
      offset = offset + mesh.numEdges();
      nodes[6] = offset + cell.vertexID(0);
      nodes[7] = offset + cell.vertexID(1);
      nodes[8] = offset + cell.vertexID(2);
      offset = offset + mesh.numVertices();
      nodes[9] = offset + cell.edgeID(0);
      nodes[10] = offset + cell.edgeID(1);
      nodes[11] = offset + cell.edgeID(2);
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
      points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
      points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
      points[6] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[7] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[8] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[9] = map(5.000000000000000e-01, 5.000000000000000e-01);
      points[10] = map(0.000000000000000e+00, 5.000000000000000e-01);
      points[11] = map(5.000000000000000e-01, 0.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
      components[4] = 0;
      components[5] = 0;
      components[6] = 1;
      components[7] = 1;
      components[8] = 1;
      components[9] = 1;
      components[10] = 1;
      components[11] = 1;
    }
  
    void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      vertex_nodes[0] = vertex;
      int offset = mesh.numVertices() + mesh.numEdges();
      vertex_nodes[1] = offset + vertex;
    }
  
    const FiniteElement& operator[] (unsigned int i) const
    {
      return *this;
    }
  
    FiniteElement& operator[] (unsigned int i)
    {
      return *this;
    }
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Vector Lagrange", "triangle", 2, 2);
      return s;
    }
    
  private:
  
    unsigned int* tensordims;
    FiniteElement** subelements;
  
  };
    
  class FunctionElement_2 : public dolfin::FiniteElement
  {
  public:
  
    FunctionElement_2() : dolfin::FiniteElement(), tensordims(0), subelements(0)
    {
      // Element is scalar, don't need to initialize tensordims
  
      // Element is simple, don't need to initialize subelements
    }
  
    ~FunctionElement_2()
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
  
    void nodemap(int nodes[], const Cell& cell, const Mesh& mesh) const
    {
      nodes[0] = cell.vertexID(0);
      nodes[1] = cell.vertexID(1);
      nodes[2] = cell.vertexID(2);
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
  
    void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      vertex_nodes[0] = vertex;
    }
  
    const FiniteElement& operator[] (unsigned int i) const
    {
      return *this;
    }
  
    FiniteElement& operator[] (unsigned int i)
    {
      return *this;
    }
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 1);
      return s;
    }
    
  private:
  
    unsigned int* tensordims;
    FiniteElement** subelements;
  
  };
  
  LinearForm(Function& w0, Function& w1, Function& w2) : dolfin::LinearForm(3)
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
    // Compute coefficients
    const real c0_0 = c[0][0];
    const real c0_1 = c[0][1];
    const real c0_2 = c[0][2];
    const real c1_0 = c[1][0];
    const real c1_1 = c[1][1];
    const real c1_2 = c[1][2];
    const real c1_3 = c[1][3];
    const real c1_4 = c[1][4];
    const real c1_5 = c[1][5];
    const real c1_6 = c[1][6];
    const real c1_7 = c[1][7];
    const real c1_8 = c[1][8];
    const real c1_9 = c[1][9];
    const real c1_10 = c[1][10];
    const real c1_11 = c[1][11];
    const real c2_0 = c[2][0];
    const real c2_1 = c[2][1];
    const real c2_2 = c[2][2];

    // Compute geometry tensors
    const real G0_0 = map.det*c0_0;
    const real G0_1 = map.det*c0_1;
    const real G0_2 = map.det*c0_2;
    const real G1_0_0_0 = map.det*c1_0*c0_0*map.g00;
    const real G1_0_0_1 = map.det*c1_0*c0_0*map.g10;
    const real G1_0_1_0 = map.det*c1_0*c0_1*map.g00;
    const real G1_0_2_1 = map.det*c1_0*c0_2*map.g10;
    const real G1_1_0_0 = map.det*c1_1*c0_0*map.g00;
    const real G1_1_0_1 = map.det*c1_1*c0_0*map.g10;
    const real G1_1_1_0 = map.det*c1_1*c0_1*map.g00;
    const real G1_1_2_1 = map.det*c1_1*c0_2*map.g10;
    const real G1_2_0_0 = map.det*c1_2*c0_0*map.g00;
    const real G1_2_0_1 = map.det*c1_2*c0_0*map.g10;
    const real G1_2_1_0 = map.det*c1_2*c0_1*map.g00;
    const real G1_2_2_1 = map.det*c1_2*c0_2*map.g10;
    const real G1_3_0_0 = map.det*c1_3*c0_0*map.g00;
    const real G1_3_0_1 = map.det*c1_3*c0_0*map.g10;
    const real G1_3_1_0 = map.det*c1_3*c0_1*map.g00;
    const real G1_3_2_1 = map.det*c1_3*c0_2*map.g10;
    const real G1_4_0_0 = map.det*c1_4*c0_0*map.g00;
    const real G1_4_0_1 = map.det*c1_4*c0_0*map.g10;
    const real G1_4_1_0 = map.det*c1_4*c0_1*map.g00;
    const real G1_4_2_1 = map.det*c1_4*c0_2*map.g10;
    const real G1_5_0_0 = map.det*c1_5*c0_0*map.g00;
    const real G1_5_0_1 = map.det*c1_5*c0_0*map.g10;
    const real G1_5_1_0 = map.det*c1_5*c0_1*map.g00;
    const real G1_5_2_1 = map.det*c1_5*c0_2*map.g10;
    const real G2_6_0_0 = map.det*c1_6*c0_0*map.g01;
    const real G2_6_0_1 = map.det*c1_6*c0_0*map.g11;
    const real G2_6_1_0 = map.det*c1_6*c0_1*map.g01;
    const real G2_6_2_1 = map.det*c1_6*c0_2*map.g11;
    const real G2_7_0_0 = map.det*c1_7*c0_0*map.g01;
    const real G2_7_0_1 = map.det*c1_7*c0_0*map.g11;
    const real G2_7_1_0 = map.det*c1_7*c0_1*map.g01;
    const real G2_7_2_1 = map.det*c1_7*c0_2*map.g11;
    const real G2_8_0_0 = map.det*c1_8*c0_0*map.g01;
    const real G2_8_0_1 = map.det*c1_8*c0_0*map.g11;
    const real G2_8_1_0 = map.det*c1_8*c0_1*map.g01;
    const real G2_8_2_1 = map.det*c1_8*c0_2*map.g11;
    const real G2_9_0_0 = map.det*c1_9*c0_0*map.g01;
    const real G2_9_0_1 = map.det*c1_9*c0_0*map.g11;
    const real G2_9_1_0 = map.det*c1_9*c0_1*map.g01;
    const real G2_9_2_1 = map.det*c1_9*c0_2*map.g11;
    const real G2_10_0_0 = map.det*c1_10*c0_0*map.g01;
    const real G2_10_0_1 = map.det*c1_10*c0_0*map.g11;
    const real G2_10_1_0 = map.det*c1_10*c0_1*map.g01;
    const real G2_10_2_1 = map.det*c1_10*c0_2*map.g11;
    const real G2_11_0_0 = map.det*c1_11*c0_0*map.g01;
    const real G2_11_0_1 = map.det*c1_11*c0_0*map.g11;
    const real G2_11_1_0 = map.det*c1_11*c0_1*map.g01;
    const real G2_11_2_1 = map.det*c1_11*c0_2*map.g11;
    const real G3_0_0_0 = map.det*c0_0*map.g00*map.g00 + map.det*c0_0*map.g01*map.g01;
    const real G3_0_0_1 = map.det*c0_0*map.g00*map.g10 + map.det*c0_0*map.g01*map.g11;
    const real G3_0_1_0 = map.det*c0_0*map.g10*map.g00 + map.det*c0_0*map.g11*map.g01;
    const real G3_0_1_1 = map.det*c0_0*map.g10*map.g10 + map.det*c0_0*map.g11*map.g11;
    const real G3_1_0_0 = map.det*c0_1*map.g00*map.g00 + map.det*c0_1*map.g01*map.g01;
    const real G3_1_1_0 = map.det*c0_1*map.g10*map.g00 + map.det*c0_1*map.g11*map.g01;
    const real G3_2_0_1 = map.det*c0_2*map.g00*map.g10 + map.det*c0_2*map.g01*map.g11;
    const real G3_2_1_1 = map.det*c0_2*map.g10*map.g10 + map.det*c0_2*map.g11*map.g11;
    const real G4_0 = map.det*c2_0;
    const real G4_1 = map.det*c2_1;
    const real G4_2 = map.det*c2_2;

    // Compute element tensor
    block[0] = 8.333333333333318e-02*G0_0 + 4.166666666666659e-02*G0_1 + 4.166666666666658e-02*G0_2 + 4.166666666666660e-04*G1_0_0_0 + 4.166666666666659e-04*G1_0_0_1 - 4.166666666666660e-04*G1_0_1_0 - 4.166666666666659e-04*G1_0_2_1 - 2.083333333333329e-04*G1_1_0_0 - 2.083333333333329e-04*G1_1_0_1 + 2.083333333333329e-04*G1_1_1_0 + 2.083333333333329e-04*G1_1_2_1 - 2.083333333333329e-04*G1_2_0_0 - 2.083333333333329e-04*G1_2_0_1 + 2.083333333333329e-04*G1_2_1_0 + 2.083333333333329e-04*G1_2_2_1 + 8.333333333333318e-04*G1_3_0_0 + 8.333333333333315e-04*G1_3_0_1 - 8.333333333333318e-04*G1_3_1_0 - 8.333333333333315e-04*G1_3_2_1 + 1.666666666666664e-03*G1_4_0_0 + 1.666666666666664e-03*G1_4_0_1 - 1.666666666666664e-03*G1_4_1_0 - 1.666666666666664e-03*G1_4_2_1 + 1.666666666666664e-03*G1_5_0_0 + 1.666666666666663e-03*G1_5_0_1 - 1.666666666666664e-03*G1_5_1_0 - 1.666666666666663e-03*G1_5_2_1 + 4.166666666666662e-04*G2_6_0_0 + 4.166666666666661e-04*G2_6_0_1 - 4.166666666666662e-04*G2_6_1_0 - 4.166666666666661e-04*G2_6_2_1 - 2.083333333333328e-04*G2_7_0_0 - 2.083333333333328e-04*G2_7_0_1 + 2.083333333333328e-04*G2_7_1_0 + 2.083333333333328e-04*G2_7_2_1 - 2.083333333333329e-04*G2_8_0_0 - 2.083333333333329e-04*G2_8_0_1 + 2.083333333333329e-04*G2_8_1_0 + 2.083333333333329e-04*G2_8_2_1 + 8.333333333333318e-04*G2_9_0_0 + 8.333333333333315e-04*G2_9_0_1 - 8.333333333333318e-04*G2_9_1_0 - 8.333333333333315e-04*G2_9_2_1 + 1.666666666666664e-03*G2_10_0_0 + 1.666666666666664e-03*G2_10_0_1 - 1.666666666666664e-03*G2_10_1_0 - 1.666666666666664e-03*G2_10_2_1 + 1.666666666666664e-03*G2_11_0_0 + 1.666666666666663e-03*G2_11_0_1 - 1.666666666666664e-03*G2_11_1_0 - 1.666666666666663e-03*G2_11_2_1 - 6.249999999999997e-05*G3_0_0_0 - 6.249999999999996e-05*G3_0_0_1 - 6.249999999999996e-05*G3_0_1_0 - 6.249999999999995e-05*G3_0_1_1 + 6.249999999999997e-05*G3_1_0_0 + 6.249999999999996e-05*G3_1_1_0 + 6.249999999999996e-05*G3_2_0_1 + 6.249999999999995e-05*G3_2_1_1 + 4.166666666666659e-03*G4_0 + 2.083333333333329e-03*G4_1 + 2.083333333333329e-03*G4_2;
    block[1] = 4.166666666666659e-02*G0_0 + 8.333333333333318e-02*G0_1 + 4.166666666666659e-02*G0_2 - 2.083333333333330e-04*G1_0_0_0 - 2.083333333333329e-04*G1_0_0_1 + 2.083333333333330e-04*G1_0_1_0 + 2.083333333333329e-04*G1_0_2_1 + 4.166666666666660e-04*G1_1_0_0 + 4.166666666666659e-04*G1_1_0_1 - 4.166666666666660e-04*G1_1_1_0 - 4.166666666666659e-04*G1_1_2_1 - 2.083333333333328e-04*G1_2_0_0 - 2.083333333333328e-04*G1_2_0_1 + 2.083333333333328e-04*G1_2_1_0 + 2.083333333333328e-04*G1_2_2_1 + 1.666666666666664e-03*G1_3_0_0 + 1.666666666666664e-03*G1_3_0_1 - 1.666666666666664e-03*G1_3_1_0 - 1.666666666666664e-03*G1_3_2_1 + 8.333333333333320e-04*G1_4_0_0 + 8.333333333333318e-04*G1_4_0_1 - 8.333333333333320e-04*G1_4_1_0 - 8.333333333333318e-04*G1_4_2_1 + 1.666666666666664e-03*G1_5_0_0 + 1.666666666666663e-03*G1_5_0_1 - 1.666666666666664e-03*G1_5_1_0 - 1.666666666666663e-03*G1_5_2_1 - 2.083333333333326e-04*G2_6_0_0 - 2.083333333333325e-04*G2_6_0_1 + 2.083333333333326e-04*G2_6_1_0 + 2.083333333333325e-04*G2_6_2_1 + 4.166666666666660e-04*G2_7_0_0 + 4.166666666666659e-04*G2_7_0_1 - 4.166666666666660e-04*G2_7_1_0 - 4.166666666666659e-04*G2_7_2_1 - 2.083333333333328e-04*G2_8_0_0 - 2.083333333333328e-04*G2_8_0_1 + 2.083333333333328e-04*G2_8_1_0 + 2.083333333333328e-04*G2_8_2_1 + 1.666666666666664e-03*G2_9_0_0 + 1.666666666666664e-03*G2_9_0_1 - 1.666666666666664e-03*G2_9_1_0 - 1.666666666666664e-03*G2_9_2_1 + 8.333333333333320e-04*G2_10_0_0 + 8.333333333333318e-04*G2_10_0_1 - 8.333333333333320e-04*G2_10_1_0 - 8.333333333333318e-04*G2_10_2_1 + 1.666666666666664e-03*G2_11_0_0 + 1.666666666666663e-03*G2_11_0_1 - 1.666666666666664e-03*G2_11_1_0 - 1.666666666666663e-03*G2_11_2_1 + 6.249999999999997e-05*G3_0_0_0 + 6.249999999999996e-05*G3_0_0_1 - 6.249999999999997e-05*G3_1_0_0 - 6.249999999999996e-05*G3_2_0_1 + 2.083333333333329e-03*G4_0 + 4.166666666666659e-03*G4_1 + 2.083333333333329e-03*G4_2;
    block[2] = 4.166666666666658e-02*G0_0 + 4.166666666666659e-02*G0_1 + 8.333333333333316e-02*G0_2 - 2.083333333333330e-04*G1_0_0_0 - 2.083333333333330e-04*G1_0_0_1 + 2.083333333333330e-04*G1_0_1_0 + 2.083333333333330e-04*G1_0_2_1 - 2.083333333333330e-04*G1_1_0_0 - 2.083333333333329e-04*G1_1_0_1 + 2.083333333333330e-04*G1_1_1_0 + 2.083333333333329e-04*G1_1_2_1 + 4.166666666666659e-04*G1_2_0_0 + 4.166666666666658e-04*G1_2_0_1 - 4.166666666666659e-04*G1_2_1_0 - 4.166666666666658e-04*G1_2_2_1 + 1.666666666666664e-03*G1_3_0_0 + 1.666666666666664e-03*G1_3_0_1 - 1.666666666666664e-03*G1_3_1_0 - 1.666666666666664e-03*G1_3_2_1 + 1.666666666666664e-03*G1_4_0_0 + 1.666666666666664e-03*G1_4_0_1 - 1.666666666666664e-03*G1_4_1_0 - 1.666666666666664e-03*G1_4_2_1 + 8.333333333333318e-04*G1_5_0_0 + 8.333333333333315e-04*G1_5_0_1 - 8.333333333333318e-04*G1_5_1_0 - 8.333333333333315e-04*G1_5_2_1 - 2.083333333333327e-04*G2_6_0_0 - 2.083333333333327e-04*G2_6_0_1 + 2.083333333333327e-04*G2_6_1_0 + 2.083333333333327e-04*G2_6_2_1 - 2.083333333333329e-04*G2_7_0_0 - 2.083333333333329e-04*G2_7_0_1 + 2.083333333333329e-04*G2_7_1_0 + 2.083333333333329e-04*G2_7_2_1 + 4.166666666666659e-04*G2_8_0_0 + 4.166666666666658e-04*G2_8_0_1 - 4.166666666666659e-04*G2_8_1_0 - 4.166666666666658e-04*G2_8_2_1 + 1.666666666666664e-03*G2_9_0_0 + 1.666666666666664e-03*G2_9_0_1 - 1.666666666666664e-03*G2_9_1_0 - 1.666666666666664e-03*G2_9_2_1 + 1.666666666666664e-03*G2_10_0_0 + 1.666666666666664e-03*G2_10_0_1 - 1.666666666666664e-03*G2_10_1_0 - 1.666666666666664e-03*G2_10_2_1 + 8.333333333333318e-04*G2_11_0_0 + 8.333333333333315e-04*G2_11_0_1 - 8.333333333333318e-04*G2_11_1_0 - 8.333333333333315e-04*G2_11_2_1 + 6.249999999999996e-05*G3_0_1_0 + 6.249999999999995e-05*G3_0_1_1 - 6.249999999999996e-05*G3_1_1_0 - 6.249999999999995e-05*G3_2_1_1 + 2.083333333333329e-03*G4_0 + 2.083333333333329e-03*G4_1 + 4.166666666666659e-03*G4_2;
  }

  // No contribution from the boundary
  void eval(real block[], const AffineMap& map, unsigned int facet) const {}   

};

} }

#endif
