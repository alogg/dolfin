// This file was automatically generated by FFC, the FEniCS Form Compiler.
// Licensed under the GNU GPL Version 2.

#ifndef __POISSONSYSTEM_BILINEAR_H
#define __POISSONSYSTEM_BILINEAR_H

#include <dolfin/NewFiniteElement.h>
#include <dolfin/BilinearForm.h>

using namespace dolfin;

/// This is the finite element for which the form is generated,
/// providing the information neccessary to do assembly.

class PoissonSystemFiniteElement : public NewFiniteElement
{
public:

  PoissonSystemFiniteElement() : NewFiniteElement() {}

  unsigned int spacedim() const
  {
    return 4;
  }

  unsigned int shapedim() const
  {
    return 3;
  }

  unsigned int tensordim(unsigned int i) const
  {
    unsigned int tensordims[] = {3};
    return tensordims[i];
  }

  unsigned int rank() const
  {
    return 1;
  }

  // FIXME: Only works for nodal basis
  unsigned int dof(unsigned int i, const Cell& cell) const
  {
    return cell.nodeID(i);
  }

  // FIXME: Only works for nodal basis
  const Point& coord(unsigned int i, const Cell& cell) const
  {
    return cell.node(i).coord();
  }

};

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class PoissonSystemBilinearForm : public BilinearForm
{
public:

  PoissonSystemBilinearForm(const NewFiniteElement& element) : BilinearForm(element) {}

  bool interior(real** A) const
  {
    // Compute geometry tensors
    real G0_00 = det*(g00*g00 + g01*g01 + g02*g02);
    real G0_01 = det*(g00*g10 + g01*g11 + g02*g12);
    real G0_02 = det*(g00*g20 + g01*g21 + g02*g22);
    real G0_10 = det*(g10*g00 + g11*g01 + g12*g02);
    real G0_11 = det*(g10*g10 + g11*g11 + g12*g12);
    real G0_12 = det*(g10*g20 + g11*g21 + g12*g22);
    real G0_20 = det*(g20*g00 + g21*g01 + g22*g02);
    real G0_21 = det*(g20*g10 + g21*g11 + g22*g12);
    real G0_22 = det*(g20*g20 + g21*g21 + g22*g22);

    // Compute element tensor
    A[0][0] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[0][1] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[0][2] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[0][3] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[0][4] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[0][5] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[0][6] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[0][7] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[0][8] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[0][9] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[0][10] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[0][11] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[1][0] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[1][1] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[1][2] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[1][3] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[1][4] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[1][5] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[1][6] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[1][7] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[1][8] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[1][9] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[1][10] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[1][11] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[2][0] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[2][1] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[2][2] = 0.166666666667*G0_00 + 0.166666666667*G0_01 + 0.166666666667*G0_02 + 0.166666666667*G0_10 + 0.166666666667*G0_11 + 0.166666666667*G0_12 + 0.166666666667*G0_20 + 0.166666666667*G0_21 + 0.166666666667*G0_22;
    A[2][3] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[2][4] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[2][5] = -0.166666666667*G0_00 - 0.166666666667*G0_10 - 0.166666666667*G0_20;
    A[2][6] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[2][7] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[2][8] = -0.166666666667*G0_01 - 0.166666666667*G0_11 - 0.166666666667*G0_21;
    A[2][9] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[2][10] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[2][11] = -0.166666666667*G0_02 - 0.166666666667*G0_12 - 0.166666666667*G0_22;
    A[3][0] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[3][1] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[3][2] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[3][3] = 0.166666666667*G0_00;
    A[3][4] = 0.166666666667*G0_00;
    A[3][5] = 0.166666666667*G0_00;
    A[3][6] = 0.166666666667*G0_01;
    A[3][7] = 0.166666666667*G0_01;
    A[3][8] = 0.166666666667*G0_01;
    A[3][9] = 0.166666666667*G0_02;
    A[3][10] = 0.166666666667*G0_02;
    A[3][11] = 0.166666666667*G0_02;
    A[4][0] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[4][1] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[4][2] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[4][3] = 0.166666666667*G0_00;
    A[4][4] = 0.166666666667*G0_00;
    A[4][5] = 0.166666666667*G0_00;
    A[4][6] = 0.166666666667*G0_01;
    A[4][7] = 0.166666666667*G0_01;
    A[4][8] = 0.166666666667*G0_01;
    A[4][9] = 0.166666666667*G0_02;
    A[4][10] = 0.166666666667*G0_02;
    A[4][11] = 0.166666666667*G0_02;
    A[5][0] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[5][1] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[5][2] = -0.166666666667*G0_00 - 0.166666666667*G0_01 - 0.166666666667*G0_02;
    A[5][3] = 0.166666666667*G0_00;
    A[5][4] = 0.166666666667*G0_00;
    A[5][5] = 0.166666666667*G0_00;
    A[5][6] = 0.166666666667*G0_01;
    A[5][7] = 0.166666666667*G0_01;
    A[5][8] = 0.166666666667*G0_01;
    A[5][9] = 0.166666666667*G0_02;
    A[5][10] = 0.166666666667*G0_02;
    A[5][11] = 0.166666666667*G0_02;
    A[6][0] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[6][1] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[6][2] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[6][3] = 0.166666666667*G0_10;
    A[6][4] = 0.166666666667*G0_10;
    A[6][5] = 0.166666666667*G0_10;
    A[6][6] = 0.166666666667*G0_11;
    A[6][7] = 0.166666666667*G0_11;
    A[6][8] = 0.166666666667*G0_11;
    A[6][9] = 0.166666666667*G0_12;
    A[6][10] = 0.166666666667*G0_12;
    A[6][11] = 0.166666666667*G0_12;
    A[7][0] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[7][1] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[7][2] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[7][3] = 0.166666666667*G0_10;
    A[7][4] = 0.166666666667*G0_10;
    A[7][5] = 0.166666666667*G0_10;
    A[7][6] = 0.166666666667*G0_11;
    A[7][7] = 0.166666666667*G0_11;
    A[7][8] = 0.166666666667*G0_11;
    A[7][9] = 0.166666666667*G0_12;
    A[7][10] = 0.166666666667*G0_12;
    A[7][11] = 0.166666666667*G0_12;
    A[8][0] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[8][1] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[8][2] = -0.166666666667*G0_10 - 0.166666666667*G0_11 - 0.166666666667*G0_12;
    A[8][3] = 0.166666666667*G0_10;
    A[8][4] = 0.166666666667*G0_10;
    A[8][5] = 0.166666666667*G0_10;
    A[8][6] = 0.166666666667*G0_11;
    A[8][7] = 0.166666666667*G0_11;
    A[8][8] = 0.166666666667*G0_11;
    A[8][9] = 0.166666666667*G0_12;
    A[8][10] = 0.166666666667*G0_12;
    A[8][11] = 0.166666666667*G0_12;
    A[9][0] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[9][1] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[9][2] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[9][3] = 0.166666666667*G0_20;
    A[9][4] = 0.166666666667*G0_20;
    A[9][5] = 0.166666666667*G0_20;
    A[9][6] = 0.166666666667*G0_21;
    A[9][7] = 0.166666666667*G0_21;
    A[9][8] = 0.166666666667*G0_21;
    A[9][9] = 0.166666666667*G0_22;
    A[9][10] = 0.166666666667*G0_22;
    A[9][11] = 0.166666666667*G0_22;
    A[10][0] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[10][1] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[10][2] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[10][3] = 0.166666666667*G0_20;
    A[10][4] = 0.166666666667*G0_20;
    A[10][5] = 0.166666666667*G0_20;
    A[10][6] = 0.166666666667*G0_21;
    A[10][7] = 0.166666666667*G0_21;
    A[10][8] = 0.166666666667*G0_21;
    A[10][9] = 0.166666666667*G0_22;
    A[10][10] = 0.166666666667*G0_22;
    A[10][11] = 0.166666666667*G0_22;
    A[11][0] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[11][1] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[11][2] = -0.166666666667*G0_20 - 0.166666666667*G0_21 - 0.166666666667*G0_22;
    A[11][3] = 0.166666666667*G0_20;
    A[11][4] = 0.166666666667*G0_20;
    A[11][5] = 0.166666666667*G0_20;
    A[11][6] = 0.166666666667*G0_21;
    A[11][7] = 0.166666666667*G0_21;
    A[11][8] = 0.166666666667*G0_21;
    A[11][9] = 0.166666666667*G0_22;
    A[11][10] = 0.166666666667*G0_22;
    A[11][11] = 0.166666666667*G0_22;

    return true;
  }

};

#endif
