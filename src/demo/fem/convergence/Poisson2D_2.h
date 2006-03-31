// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.1.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __POISSON2D_2_H
#define __POISSON2D_2_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/Vector.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/FiniteElementSpec.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace Poisson2D_2 {

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
      return 6;
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
      int offset = mesh.numVertices();
      nodes[3] = offset + cell.edgeID(0);
      nodes[4] = offset + cell.edgeID(1);
      nodes[5] = offset + cell.edgeID(2);
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
      points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
      points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
      components[4] = 0;
      components[5] = 0;
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
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 2);
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
      return 6;
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
      int offset = mesh.numVertices();
      nodes[3] = offset + cell.edgeID(0);
      nodes[4] = offset + cell.edgeID(1);
      nodes[5] = offset + cell.edgeID(2);
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
      points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
      points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
      components[4] = 0;
      components[5] = 0;
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
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 2);
      return s;
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
    const real G0_0_0 = map.det*(map.g00*map.g00 + map.g01*map.g01);
    const real G0_0_1 = map.det*(map.g00*map.g10 + map.g01*map.g11);
    const real G0_1_0 = map.det*(map.g10*map.g00 + map.g11*map.g01);
    const real G0_1_1 = map.det*(map.g10*map.g10 + map.g11*map.g11);

    // Compute element tensor
    block[9] = 6.666666666666654e-01*G0_0_1;
    block[18] = block[9] + -6.666666666666652e-01*G0_0_1;
    block[34] = block[9] + 6.666666666666654e-01*G0_1_0;
    block[13] = block[18] + -1.666666666666663e-01*G0_1_0;
    block[3] = block[18];
    block[17] = block[18];
    block[8] = block[9] + -8.333333333333318e-01*G0_0_1;
    block[11] = -block[9] + -6.666666666666654e-01*G0_0_0;
    block[30] = block[11];
    block[14] = block[18] + 4.999999999999991e-01*G0_1_1;
    block[25] = block[18];
    block[7] = block[18] + 4.999999999999989e-01*G0_0_0;
    block[20] = block[9];
    block[0] = block[14] + 4.999999999999992e-01*G0_0_0 + 4.999999999999991e-01*G0_0_1 + 4.999999999999991e-01*G0_1_0;
    block[29] = block[34];
    block[4] = -block[9] + -6.666666666666656e-01*G0_1_1;
    block[10] = block[18];
    block[26] = block[4];
    block[6] = -block[8] + 1.666666666666663e-01*G0_0_0;
    block[23] = -block[34] + -1.333333333333331e+00*G0_1_1;
    block[21] = -block[23] + 1.333333333333330e+00*G0_0_0;
    block[28] = block[21];
    block[33] = block[23];
    block[35] = block[21];
    block[1] = -block[13] + 1.666666666666663e-01*G0_0_0;
    block[31] = -4.000000000000001e+00*block[1];
    block[5] = block[31];
    block[32] = block[18];
    block[19] = block[18] + 6.666666666666655e-01*G0_1_0;
    block[24] = -block[19] + -6.666666666666656e-01*G0_1_1;
    block[2] = -block[8] + 1.666666666666665e-01*G0_1_1;
    block[12] = -block[13] + 1.666666666666665e-01*G0_1_1;
    block[27] = -block[34] + -1.333333333333330e+00*G0_0_0;
    block[15] = block[19];
    block[22] = block[27];
    block[16] = block[24];
  }

  // No contribution from the boundary
  void eval(real block[], const AffineMap& map, unsigned int boundary) const {}   

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
      return 6;
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
      int offset = mesh.numVertices();
      nodes[3] = offset + cell.edgeID(0);
      nodes[4] = offset + cell.edgeID(1);
      nodes[5] = offset + cell.edgeID(2);
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
      points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
      points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
      components[4] = 0;
      components[5] = 0;
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
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 2);
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
      return 6;
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
      int offset = mesh.numVertices();
      nodes[3] = offset + cell.edgeID(0);
      nodes[4] = offset + cell.edgeID(1);
      nodes[5] = offset + cell.edgeID(2);
    }
  
    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
      points[3] = map(5.000000000000000e-01, 5.000000000000000e-01);
      points[4] = map(0.000000000000000e+00, 5.000000000000000e-01);
      points[5] = map(5.000000000000000e-01, 0.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
      components[4] = 0;
      components[5] = 0;
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
  
    FiniteElementSpec spec() const
    {
      FiniteElementSpec s("Lagrange", "triangle", 2);
      return s;
    }
    
  private:
  
    unsigned int* tensordims;
    FiniteElement** subelements;
  
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
    // Compute coefficients
    const real c0_0 = c[0][0];
    const real c0_1 = c[0][1];
    const real c0_2 = c[0][2];
    const real c0_3 = c[0][3];
    const real c0_4 = c[0][4];
    const real c0_5 = c[0][5];

    // Compute geometry tensors
    const real G0_0 = map.det*c0_0;
    const real G0_1 = map.det*c0_1;
    const real G0_2 = map.det*c0_2;
    const real G0_3 = map.det*c0_3;
    const real G0_4 = map.det*c0_4;
    const real G0_5 = map.det*c0_5;

    // Compute element tensor
    block[0] = 1.666666666666665e-02*G0_0 - 2.777777777777774e-03*G0_1 - 2.777777777777775e-03*G0_2 - 1.111111111111110e-02*G0_3;
    block[1] = -2.777777777777774e-03*G0_0 + 1.666666666666665e-02*G0_1 - 2.777777777777776e-03*G0_2 - 1.111111111111111e-02*G0_4;
    block[2] = -2.777777777777775e-03*G0_0 - 2.777777777777776e-03*G0_1 + 1.666666666666666e-02*G0_2 - 1.111111111111111e-02*G0_5;
    block[3] = -1.111111111111110e-02*G0_0 + 8.888888888888882e-02*G0_3 + 4.444444444444443e-02*G0_4 + 4.444444444444443e-02*G0_5;
    block[4] = -1.111111111111111e-02*G0_1 + 4.444444444444443e-02*G0_3 + 8.888888888888884e-02*G0_4 + 4.444444444444442e-02*G0_5;
    block[5] = -1.111111111111111e-02*G0_2 + 4.444444444444443e-02*G0_3 + 4.444444444444443e-02*G0_4 + 8.888888888888882e-02*G0_5;
  }

  // No contribution from the boundary
  void eval(real block[], const AffineMap& map, unsigned int boundary) const {}   

};

} }

#endif
