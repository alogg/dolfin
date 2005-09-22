// Automatically generated by FFC, the FEniCS Form Compiler, version 0.1.9.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __ELASTICITYUPDATEDMASS_H
#define __ELASTICITYUPDATEDMASS_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/Vector.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace ElasticityUpdatedMass {

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

    void vertexeval(real values[], unsigned int vertex, const Vector& x, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      values[0] = x(vertex);
      int offset = mesh.noNodes();
      values[1] = x(offset + vertex);
      offset = offset + mesh.noNodes();
      values[2] = x(offset + vertex);
    }

  private:

    unsigned int* tensordims;

  };
    
  class TrialElement : public dolfin::FiniteElement
  {
  public:

    TrialElement() : dolfin::FiniteElement(), tensordims(0)
    {
      tensordims = new unsigned int [1];
      tensordims[0] = 3;
    }

    ~TrialElement()
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

    void vertexeval(real values[], unsigned int vertex, const Vector& x, const Mesh& mesh) const
    {
      // FIXME: Temporary fix for Lagrange elements
      values[0] = x(vertex);
      int offset = mesh.noNodes();
      values[1] = x(offset + vertex);
      offset = offset + mesh.noNodes();
      values[2] = x(offset + vertex);
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
    // Compute geometry tensors
    real G0_0 = map.det*c[0][0];
    real G0_1 = map.det*c[0][1];
    real G0_2 = map.det*c[0][2];
    real G0_3 = map.det*c[0][3];

    // Compute element tensor
    block[0] = 8.333333333333307e-03*G0_0 + 2.777777777777771e-03*G0_1 + 2.777777777777770e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[1] = 2.777777777777769e-03*G0_0 + 2.777777777777770e-03*G0_1 + 1.388888888888885e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[2] = 2.777777777777769e-03*G0_0 + 1.388888888888885e-03*G0_1 + 2.777777777777770e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[3] = 2.777777777777770e-03*G0_0 + 1.388888888888885e-03*G0_1 + 1.388888888888886e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[4] = 0.000000000000000e+00;
    block[5] = 0.000000000000000e+00;
    block[6] = 0.000000000000000e+00;
    block[7] = 0.000000000000000e+00;
    block[8] = 0.000000000000000e+00;
    block[9] = 0.000000000000000e+00;
    block[10] = 0.000000000000000e+00;
    block[11] = 0.000000000000000e+00;
    block[12] = 2.777777777777769e-03*G0_0 + 2.777777777777770e-03*G0_1 + 1.388888888888885e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[13] = 2.777777777777769e-03*G0_0 + 8.333333333333307e-03*G0_1 + 2.777777777777770e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[14] = 1.388888888888884e-03*G0_0 + 2.777777777777769e-03*G0_1 + 2.777777777777769e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[15] = 1.388888888888885e-03*G0_0 + 2.777777777777770e-03*G0_1 + 1.388888888888885e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[16] = 0.000000000000000e+00;
    block[17] = 0.000000000000000e+00;
    block[18] = 0.000000000000000e+00;
    block[19] = 0.000000000000000e+00;
    block[20] = 0.000000000000000e+00;
    block[21] = 0.000000000000000e+00;
    block[22] = 0.000000000000000e+00;
    block[23] = 0.000000000000000e+00;
    block[24] = 2.777777777777770e-03*G0_0 + 1.388888888888885e-03*G0_1 + 2.777777777777770e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[25] = 1.388888888888884e-03*G0_0 + 2.777777777777770e-03*G0_1 + 2.777777777777769e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[26] = 2.777777777777768e-03*G0_0 + 2.777777777777770e-03*G0_1 + 8.333333333333309e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[27] = 1.388888888888885e-03*G0_0 + 1.388888888888885e-03*G0_1 + 2.777777777777771e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[28] = 0.000000000000000e+00;
    block[29] = 0.000000000000000e+00;
    block[30] = 0.000000000000000e+00;
    block[31] = 0.000000000000000e+00;
    block[32] = 0.000000000000000e+00;
    block[33] = 0.000000000000000e+00;
    block[34] = 0.000000000000000e+00;
    block[35] = 0.000000000000000e+00;
    block[36] = 2.777777777777770e-03*G0_0 + 1.388888888888885e-03*G0_1 + 1.388888888888886e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[37] = 1.388888888888885e-03*G0_0 + 2.777777777777770e-03*G0_1 + 1.388888888888885e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[38] = 1.388888888888884e-03*G0_0 + 1.388888888888885e-03*G0_1 + 2.777777777777771e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[39] = 2.777777777777770e-03*G0_0 + 2.777777777777771e-03*G0_1 + 2.777777777777771e-03*G0_2 + 8.333333333333311e-03*G0_3;
    block[40] = 0.000000000000000e+00;
    block[41] = 0.000000000000000e+00;
    block[42] = 0.000000000000000e+00;
    block[43] = 0.000000000000000e+00;
    block[44] = 0.000000000000000e+00;
    block[45] = 0.000000000000000e+00;
    block[46] = 0.000000000000000e+00;
    block[47] = 0.000000000000000e+00;
    block[48] = 0.000000000000000e+00;
    block[49] = 0.000000000000000e+00;
    block[50] = 0.000000000000000e+00;
    block[51] = 0.000000000000000e+00;
    block[52] = 8.333333333333304e-03*G0_0 + 2.777777777777769e-03*G0_1 + 2.777777777777769e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[53] = 2.777777777777768e-03*G0_0 + 2.777777777777769e-03*G0_1 + 1.388888888888884e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[54] = 2.777777777777770e-03*G0_0 + 1.388888888888885e-03*G0_1 + 2.777777777777769e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[55] = 2.777777777777769e-03*G0_0 + 1.388888888888885e-03*G0_1 + 1.388888888888885e-03*G0_2 + 2.777777777777769e-03*G0_3;
    block[56] = 0.000000000000000e+00;
    block[57] = 0.000000000000000e+00;
    block[58] = 0.000000000000000e+00;
    block[59] = 0.000000000000000e+00;
    block[60] = 0.000000000000000e+00;
    block[61] = 0.000000000000000e+00;
    block[62] = 0.000000000000000e+00;
    block[63] = 0.000000000000000e+00;
    block[64] = 2.777777777777768e-03*G0_0 + 2.777777777777768e-03*G0_1 + 1.388888888888884e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[65] = 2.777777777777769e-03*G0_0 + 8.333333333333307e-03*G0_1 + 2.777777777777771e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[66] = 1.388888888888885e-03*G0_0 + 2.777777777777771e-03*G0_1 + 2.777777777777771e-03*G0_2 + 1.388888888888886e-03*G0_3;
    block[67] = 1.388888888888885e-03*G0_0 + 2.777777777777770e-03*G0_1 + 1.388888888888886e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[68] = 0.000000000000000e+00;
    block[69] = 0.000000000000000e+00;
    block[70] = 0.000000000000000e+00;
    block[71] = 0.000000000000000e+00;
    block[72] = 0.000000000000000e+00;
    block[73] = 0.000000000000000e+00;
    block[74] = 0.000000000000000e+00;
    block[75] = 0.000000000000000e+00;
    block[76] = 2.777777777777770e-03*G0_0 + 1.388888888888885e-03*G0_1 + 2.777777777777769e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[77] = 1.388888888888885e-03*G0_0 + 2.777777777777771e-03*G0_1 + 2.777777777777771e-03*G0_2 + 1.388888888888886e-03*G0_3;
    block[78] = 2.777777777777770e-03*G0_0 + 2.777777777777771e-03*G0_1 + 8.333333333333311e-03*G0_2 + 2.777777777777772e-03*G0_3;
    block[79] = 1.388888888888885e-03*G0_0 + 1.388888888888886e-03*G0_1 + 2.777777777777771e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[80] = 0.000000000000000e+00;
    block[81] = 0.000000000000000e+00;
    block[82] = 0.000000000000000e+00;
    block[83] = 0.000000000000000e+00;
    block[84] = 0.000000000000000e+00;
    block[85] = 0.000000000000000e+00;
    block[86] = 0.000000000000000e+00;
    block[87] = 0.000000000000000e+00;
    block[88] = 2.777777777777769e-03*G0_0 + 1.388888888888885e-03*G0_1 + 1.388888888888885e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[89] = 1.388888888888885e-03*G0_0 + 2.777777777777770e-03*G0_1 + 1.388888888888885e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[90] = 1.388888888888885e-03*G0_0 + 1.388888888888886e-03*G0_1 + 2.777777777777771e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[91] = 2.777777777777770e-03*G0_0 + 2.777777777777771e-03*G0_1 + 2.777777777777771e-03*G0_2 + 8.333333333333311e-03*G0_3;
    block[92] = 0.000000000000000e+00;
    block[93] = 0.000000000000000e+00;
    block[94] = 0.000000000000000e+00;
    block[95] = 0.000000000000000e+00;
    block[96] = 0.000000000000000e+00;
    block[97] = 0.000000000000000e+00;
    block[98] = 0.000000000000000e+00;
    block[99] = 0.000000000000000e+00;
    block[100] = 0.000000000000000e+00;
    block[101] = 0.000000000000000e+00;
    block[102] = 0.000000000000000e+00;
    block[103] = 0.000000000000000e+00;
    block[104] = 8.333333333333304e-03*G0_0 + 2.777777777777769e-03*G0_1 + 2.777777777777769e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[105] = 2.777777777777768e-03*G0_0 + 2.777777777777769e-03*G0_1 + 1.388888888888884e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[106] = 2.777777777777770e-03*G0_0 + 1.388888888888885e-03*G0_1 + 2.777777777777769e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[107] = 2.777777777777769e-03*G0_0 + 1.388888888888885e-03*G0_1 + 1.388888888888885e-03*G0_2 + 2.777777777777769e-03*G0_3;
    block[108] = 0.000000000000000e+00;
    block[109] = 0.000000000000000e+00;
    block[110] = 0.000000000000000e+00;
    block[111] = 0.000000000000000e+00;
    block[112] = 0.000000000000000e+00;
    block[113] = 0.000000000000000e+00;
    block[114] = 0.000000000000000e+00;
    block[115] = 0.000000000000000e+00;
    block[116] = 2.777777777777768e-03*G0_0 + 2.777777777777768e-03*G0_1 + 1.388888888888884e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[117] = 2.777777777777769e-03*G0_0 + 8.333333333333307e-03*G0_1 + 2.777777777777771e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[118] = 1.388888888888885e-03*G0_0 + 2.777777777777771e-03*G0_1 + 2.777777777777771e-03*G0_2 + 1.388888888888886e-03*G0_3;
    block[119] = 1.388888888888885e-03*G0_0 + 2.777777777777770e-03*G0_1 + 1.388888888888886e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[120] = 0.000000000000000e+00;
    block[121] = 0.000000000000000e+00;
    block[122] = 0.000000000000000e+00;
    block[123] = 0.000000000000000e+00;
    block[124] = 0.000000000000000e+00;
    block[125] = 0.000000000000000e+00;
    block[126] = 0.000000000000000e+00;
    block[127] = 0.000000000000000e+00;
    block[128] = 2.777777777777770e-03*G0_0 + 1.388888888888885e-03*G0_1 + 2.777777777777769e-03*G0_2 + 1.388888888888885e-03*G0_3;
    block[129] = 1.388888888888885e-03*G0_0 + 2.777777777777771e-03*G0_1 + 2.777777777777771e-03*G0_2 + 1.388888888888886e-03*G0_3;
    block[130] = 2.777777777777770e-03*G0_0 + 2.777777777777771e-03*G0_1 + 8.333333333333311e-03*G0_2 + 2.777777777777772e-03*G0_3;
    block[131] = 1.388888888888885e-03*G0_0 + 1.388888888888886e-03*G0_1 + 2.777777777777771e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[132] = 0.000000000000000e+00;
    block[133] = 0.000000000000000e+00;
    block[134] = 0.000000000000000e+00;
    block[135] = 0.000000000000000e+00;
    block[136] = 0.000000000000000e+00;
    block[137] = 0.000000000000000e+00;
    block[138] = 0.000000000000000e+00;
    block[139] = 0.000000000000000e+00;
    block[140] = 2.777777777777769e-03*G0_0 + 1.388888888888885e-03*G0_1 + 1.388888888888885e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[141] = 1.388888888888885e-03*G0_0 + 2.777777777777770e-03*G0_1 + 1.388888888888885e-03*G0_2 + 2.777777777777770e-03*G0_3;
    block[142] = 1.388888888888885e-03*G0_0 + 1.388888888888886e-03*G0_1 + 2.777777777777771e-03*G0_2 + 2.777777777777771e-03*G0_3;
    block[143] = 2.777777777777770e-03*G0_0 + 2.777777777777771e-03*G0_1 + 2.777777777777771e-03*G0_2 + 8.333333333333311e-03*G0_3;
  }

};

} }

#endif
