// Automatically generated by FFC, the FEniCS Form Compiler, version 0.1.7.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __./SRC/MODULES/ELASTICITY/DOLFIN/ELASTICITYMASS_H
#define __./SRC/MODULES/ELASTICITY/DOLFIN/ELASTICITYMASS_H

#include <dolfin/FiniteElement.h>
#include <dolfin/LinearForm.h>
#include <dolfin/BilinearForm.h>

namespace dolfin { namespace ./src/modules/elasticity/dolfin/ElasticityMass {

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

    // FIXME: Only works for nodal basis
    inline unsigned int dof(unsigned int i, const Cell& cell, const Mesh& mesh) const
    {
      return (i/4) * mesh.noNodes() + cell.nodeID(i % 4);
    }

    // FIXME: Only works for nodal basis
    inline const Point coord(unsigned int i, const Cell& cell, const Mesh& mesh) const
    {
      return cell.node(i % 4).coord();
    }

    // FIXME: New version replacing dof()
    inline void dofmap(int dofs[], const Cell& cell, const Mesh& mesh) const
    {
      dofs[0] = cell.nodeID(0);
      dofs[1] = cell.nodeID(1);
      dofs[2] = cell.nodeID(2);
      dofs[3] = cell.nodeID(3);
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

    // FIXME: Only works for nodal basis
    inline unsigned int dof(unsigned int i, const Cell& cell, const Mesh& mesh) const
    {
      return (i/4) * mesh.noNodes() + cell.nodeID(i % 4);
    }

    // FIXME: Only works for nodal basis
    inline const Point coord(unsigned int i, const Cell& cell, const Mesh& mesh) const
    {
      return cell.node(i % 4).coord();
    }

    // FIXME: New version replacing dof()
    inline void dofmap(int dofs[], const Cell& cell, const Mesh& mesh) const
    {
      dofs[0] = cell.nodeID(0);
      dofs[1] = cell.nodeID(1);
      dofs[2] = cell.nodeID(2);
      dofs[3] = cell.nodeID(3);
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

  bool interior(real* block) const
  {
    // Compute geometry tensors
    real G0_ = det;

    // Compute element tensor
    block[0] = 1.666666666666662e-02*G0_;
    block[1] = 8.333333333333309e-03*G0_;
    block[2] = 8.333333333333309e-03*G0_;
    block[3] = 8.333333333333311e-03*G0_;
    block[4] = 0.000000000000000e+00;
    block[5] = 0.000000000000000e+00;
    block[6] = 0.000000000000000e+00;
    block[7] = 0.000000000000000e+00;
    block[8] = 0.000000000000000e+00;
    block[9] = 0.000000000000000e+00;
    block[10] = 0.000000000000000e+00;
    block[11] = 0.000000000000000e+00;
    block[12] = 8.333333333333309e-03*G0_;
    block[13] = 1.666666666666661e-02*G0_;
    block[14] = 8.333333333333307e-03*G0_;
    block[15] = 8.333333333333307e-03*G0_;
    block[16] = 0.000000000000000e+00;
    block[17] = 0.000000000000000e+00;
    block[18] = 0.000000000000000e+00;
    block[19] = 0.000000000000000e+00;
    block[20] = 0.000000000000000e+00;
    block[21] = 0.000000000000000e+00;
    block[22] = 0.000000000000000e+00;
    block[23] = 0.000000000000000e+00;
    block[24] = 8.333333333333307e-03*G0_;
    block[25] = 8.333333333333309e-03*G0_;
    block[26] = 1.666666666666661e-02*G0_;
    block[27] = 8.333333333333311e-03*G0_;
    block[28] = 0.000000000000000e+00;
    block[29] = 0.000000000000000e+00;
    block[30] = 0.000000000000000e+00;
    block[31] = 0.000000000000000e+00;
    block[32] = 0.000000000000000e+00;
    block[33] = 0.000000000000000e+00;
    block[34] = 0.000000000000000e+00;
    block[35] = 0.000000000000000e+00;
    block[36] = 8.333333333333311e-03*G0_;
    block[37] = 8.333333333333307e-03*G0_;
    block[38] = 8.333333333333311e-03*G0_;
    block[39] = 1.666666666666662e-02*G0_;
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
    block[52] = 1.666666666666661e-02*G0_;
    block[53] = 8.333333333333307e-03*G0_;
    block[54] = 8.333333333333307e-03*G0_;
    block[55] = 8.333333333333305e-03*G0_;
    block[56] = 0.000000000000000e+00;
    block[57] = 0.000000000000000e+00;
    block[58] = 0.000000000000000e+00;
    block[59] = 0.000000000000000e+00;
    block[60] = 0.000000000000000e+00;
    block[61] = 0.000000000000000e+00;
    block[62] = 0.000000000000000e+00;
    block[63] = 0.000000000000000e+00;
    block[64] = 8.333333333333307e-03*G0_;
    block[65] = 1.666666666666662e-02*G0_;
    block[66] = 8.333333333333311e-03*G0_;
    block[67] = 8.333333333333311e-03*G0_;
    block[68] = 0.000000000000000e+00;
    block[69] = 0.000000000000000e+00;
    block[70] = 0.000000000000000e+00;
    block[71] = 0.000000000000000e+00;
    block[72] = 0.000000000000000e+00;
    block[73] = 0.000000000000000e+00;
    block[74] = 0.000000000000000e+00;
    block[75] = 0.000000000000000e+00;
    block[76] = 8.333333333333307e-03*G0_;
    block[77] = 8.333333333333311e-03*G0_;
    block[78] = 1.666666666666662e-02*G0_;
    block[79] = 8.333333333333314e-03*G0_;
    block[80] = 0.000000000000000e+00;
    block[81] = 0.000000000000000e+00;
    block[82] = 0.000000000000000e+00;
    block[83] = 0.000000000000000e+00;
    block[84] = 0.000000000000000e+00;
    block[85] = 0.000000000000000e+00;
    block[86] = 0.000000000000000e+00;
    block[87] = 0.000000000000000e+00;
    block[88] = 8.333333333333305e-03*G0_;
    block[89] = 8.333333333333311e-03*G0_;
    block[90] = 8.333333333333314e-03*G0_;
    block[91] = 1.666666666666662e-02*G0_;
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
    block[104] = 1.666666666666661e-02*G0_;
    block[105] = 8.333333333333307e-03*G0_;
    block[106] = 8.333333333333307e-03*G0_;
    block[107] = 8.333333333333305e-03*G0_;
    block[108] = 0.000000000000000e+00;
    block[109] = 0.000000000000000e+00;
    block[110] = 0.000000000000000e+00;
    block[111] = 0.000000000000000e+00;
    block[112] = 0.000000000000000e+00;
    block[113] = 0.000000000000000e+00;
    block[114] = 0.000000000000000e+00;
    block[115] = 0.000000000000000e+00;
    block[116] = 8.333333333333307e-03*G0_;
    block[117] = 1.666666666666662e-02*G0_;
    block[118] = 8.333333333333311e-03*G0_;
    block[119] = 8.333333333333311e-03*G0_;
    block[120] = 0.000000000000000e+00;
    block[121] = 0.000000000000000e+00;
    block[122] = 0.000000000000000e+00;
    block[123] = 0.000000000000000e+00;
    block[124] = 0.000000000000000e+00;
    block[125] = 0.000000000000000e+00;
    block[126] = 0.000000000000000e+00;
    block[127] = 0.000000000000000e+00;
    block[128] = 8.333333333333307e-03*G0_;
    block[129] = 8.333333333333311e-03*G0_;
    block[130] = 1.666666666666662e-02*G0_;
    block[131] = 8.333333333333312e-03*G0_;
    block[132] = 0.000000000000000e+00;
    block[133] = 0.000000000000000e+00;
    block[134] = 0.000000000000000e+00;
    block[135] = 0.000000000000000e+00;
    block[136] = 0.000000000000000e+00;
    block[137] = 0.000000000000000e+00;
    block[138] = 0.000000000000000e+00;
    block[139] = 0.000000000000000e+00;
    block[140] = 8.333333333333305e-03*G0_;
    block[141] = 8.333333333333311e-03*G0_;
    block[142] = 8.333333333333312e-03*G0_;
    block[143] = 1.666666666666662e-02*G0_;

    return true;
  }

};

} }

#endif
