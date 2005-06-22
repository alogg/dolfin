// Automatically generated by FFC, the FEniCS Form Compiler, version 0.1.8.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __STIFFNESSMATRIX3D_H
#define __STIFFNESSMATRIX3D_H

#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace StiffnessMatrix3D {

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

  private:

    unsigned int* tensordims;

  };

  BilinearForm(const real& c0) : dolfin::BilinearForm(0), c0(c0)
  {
    // Create finite element for test space
    _test = new TestElement();

    // Create finite element for trial space
    _trial = new TrialElement();
  }

  void eval(real block[], const AffineMap& map) const
  {
    // Compute geometry tensors
    real G0_0_0 = map.det*c0*(map.g00*map.g00 + map.g01*map.g01 + map.g02*map.g02);
    real G0_0_1 = map.det*c0*(map.g00*map.g10 + map.g01*map.g11 + map.g02*map.g12);
    real G0_0_2 = map.det*c0*(map.g00*map.g20 + map.g01*map.g21 + map.g02*map.g22);
    real G0_1_0 = map.det*c0*(map.g10*map.g00 + map.g11*map.g01 + map.g12*map.g02);
    real G0_1_1 = map.det*c0*(map.g10*map.g10 + map.g11*map.g11 + map.g12*map.g12);
    real G0_1_2 = map.det*c0*(map.g10*map.g20 + map.g11*map.g21 + map.g12*map.g22);
    real G0_2_0 = map.det*c0*(map.g20*map.g00 + map.g21*map.g01 + map.g22*map.g02);
    real G0_2_1 = map.det*c0*(map.g20*map.g10 + map.g21*map.g11 + map.g22*map.g12);
    real G0_2_2 = map.det*c0*(map.g20*map.g20 + map.g21*map.g21 + map.g22*map.g22);

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
        
private:

  const real& c0;

};

} }

#endif
