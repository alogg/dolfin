// Automatically generated by FFC, the FEniCS Form Compiler, version 0.1.9.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __POISSON3D_1_H
#define __POISSON3D_1_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/Vector.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace Poisson3D_1 {

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
      return 4;
    }

    inline unsigned int shapedim() const
    {
      return 3;
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
      dofs[3] = cell.nodeID(3);
    }

    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
      points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
    }

    void vertexeval(real values[], unsigned int vertex, const Vector& x, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      values[0] = x(vertex);
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
      return 4;
    }

    inline unsigned int shapedim() const
    {
      return 3;
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
      dofs[3] = cell.nodeID(3);
    }

    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
      points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
    }

    void vertexeval(real values[], unsigned int vertex, const Vector& x, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      values[0] = x(vertex);
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
    real G0_0_0 = map.det*(map.g00*map.g00 + map.g01*map.g01 + map.g02*map.g02);
    real G0_0_1 = map.det*(map.g00*map.g10 + map.g01*map.g11 + map.g02*map.g12);
    real G0_0_2 = map.det*(map.g00*map.g20 + map.g01*map.g21 + map.g02*map.g22);
    real G0_1_0 = map.det*(map.g10*map.g00 + map.g11*map.g01 + map.g12*map.g02);
    real G0_1_1 = map.det*(map.g10*map.g10 + map.g11*map.g11 + map.g12*map.g12);
    real G0_1_2 = map.det*(map.g10*map.g20 + map.g11*map.g21 + map.g12*map.g22);
    real G0_2_0 = map.det*(map.g20*map.g00 + map.g21*map.g01 + map.g22*map.g02);
    real G0_2_1 = map.det*(map.g20*map.g10 + map.g21*map.g11 + map.g22*map.g12);
    real G0_2_2 = map.det*(map.g20*map.g20 + map.g21*map.g21 + map.g22*map.g22);

    // Compute element tensor
    block[0] = 1.666666666666665e-01*G0_0_0 + 1.666666666666665e-01*G0_0_1 + 1.666666666666665e-01*G0_0_2 + 1.666666666666665e-01*G0_1_0 + 1.666666666666665e-01*G0_1_1 + 1.666666666666665e-01*G0_1_2 + 1.666666666666665e-01*G0_2_0 + 1.666666666666665e-01*G0_2_1 + 1.666666666666665e-01*G0_2_2;
    block[1] = -1.666666666666665e-01*G0_0_0 - 1.666666666666665e-01*G0_1_0 - 1.666666666666665e-01*G0_2_0;
    block[2] = -1.666666666666665e-01*G0_0_1 - 1.666666666666665e-01*G0_1_1 - 1.666666666666665e-01*G0_2_1;
    block[3] = -1.666666666666665e-01*G0_0_2 - 1.666666666666665e-01*G0_1_2 - 1.666666666666665e-01*G0_2_2;
    block[4] = -1.666666666666665e-01*G0_0_0 - 1.666666666666665e-01*G0_0_1 - 1.666666666666665e-01*G0_0_2;
    block[5] = 1.666666666666665e-01*G0_0_0;
    block[6] = 1.666666666666665e-01*G0_0_1;
    block[7] = 1.666666666666665e-01*G0_0_2;
    block[8] = -1.666666666666665e-01*G0_1_0 - 1.666666666666665e-01*G0_1_1 - 1.666666666666665e-01*G0_1_2;
    block[9] = 1.666666666666665e-01*G0_1_0;
    block[10] = 1.666666666666665e-01*G0_1_1;
    block[11] = 1.666666666666665e-01*G0_1_2;
    block[12] = -1.666666666666665e-01*G0_2_0 - 1.666666666666665e-01*G0_2_1 - 1.666666666666665e-01*G0_2_2;
    block[13] = 1.666666666666665e-01*G0_2_0;
    block[14] = 1.666666666666665e-01*G0_2_1;
    block[15] = 1.666666666666665e-01*G0_2_2;
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
      return 4;
    }

    inline unsigned int shapedim() const
    {
      return 3;
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
      dofs[3] = cell.nodeID(3);
    }

    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
      points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
    }

    void vertexeval(real values[], unsigned int vertex, const Vector& x, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      values[0] = x(vertex);
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
      return 4;
    }

    inline unsigned int shapedim() const
    {
      return 3;
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
      dofs[3] = cell.nodeID(3);
    }

    void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
    {
      points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
      points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
      points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
      components[0] = 0;
      components[1] = 0;
      components[2] = 0;
      components[3] = 0;
    }

    void vertexeval(real values[], unsigned int vertex, const Vector& x, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      values[0] = x(vertex);
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
    real G0_3 = map.det*c[0][3];

    // Compute element tensor
    block[0] = 1.666666666666662e-02*G0_0 + 8.333333333333309e-03*G0_1 + 8.333333333333311e-03*G0_2 + 8.333333333333312e-03*G0_3;
    block[1] = 8.333333333333309e-03*G0_0 + 1.666666666666662e-02*G0_1 + 8.333333333333312e-03*G0_2 + 8.333333333333312e-03*G0_3;
    block[2] = 8.333333333333307e-03*G0_0 + 8.333333333333311e-03*G0_1 + 1.666666666666662e-02*G0_2 + 8.333333333333314e-03*G0_3;
    block[3] = 8.333333333333312e-03*G0_0 + 8.333333333333312e-03*G0_1 + 8.333333333333314e-03*G0_2 + 1.666666666666662e-02*G0_3;
  }

};

} }

#endif
