// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.5.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __POISSON2D_3_H
#define __POISSON2D_3_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/FiniteElementSpec.h>
#include <dolfin/BilinearForm.h>
#include <dolfin/LinearForm.h>
#include <dolfin/Functional.h>
#include <dolfin/FEM.h>

namespace dolfin { namespace Poisson2D_3 {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
{
public:

  class TestElement;

  class TrialElement;

  BilinearForm();
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map, real det) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, real det, unsigned int facet) const;

  bool interior_boundary_contribution() const;

  void eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const;

};

class BilinearForm::TestElement : public dolfin::FiniteElement
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
    return 10;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    error("Element is scalar.");
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
    static unsigned int edge_reordering_0[2][2] = {{0, 1}, {1, 0}};
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    int alignment = cell.alignment(1, 0);
    int offset = mesh.topology().size(0);
    nodes[3] = offset + 2*cell.entities(1)[0] + edge_reordering_0[alignment][0];
    nodes[4] = offset + 2*cell.entities(1)[0] + edge_reordering_0[alignment][1];
    alignment = cell.alignment(1, 1);
    nodes[5] = offset + 2*cell.entities(1)[1] + edge_reordering_0[alignment][0];
    nodes[6] = offset + 2*cell.entities(1)[1] + edge_reordering_0[alignment][1];
    alignment = cell.alignment(1, 2);
    nodes[7] = offset + 2*cell.entities(1)[2] + edge_reordering_0[alignment][0];
    nodes[8] = offset + 2*cell.entities(1)[2] + edge_reordering_0[alignment][1];
    offset = offset + 2*mesh.topology().size(1);
    nodes[9] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(6.666666666666667e-01, 3.333333333333333e-01);
    points[4] = map(3.333333333333334e-01, 6.666666666666666e-01);
    points[5] = map(0.000000000000000e+00, 6.666666666666667e-01);
    points[6] = map(0.000000000000000e+00, 3.333333333333334e-01);
    points[7] = map(3.333333333333333e-01, 0.000000000000000e+00);
    points[8] = map(6.666666666666666e-01, 0.000000000000000e+00);
    points[9] = map(3.333333333333333e-01, 3.333333333333333e-01);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 0;
    components[5] = 0;
    components[6] = 0;
    components[7] = 0;
    components[8] = 0;
    components[9] = 0;
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
    FiniteElementSpec s("Lagrange", "triangle", 3);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class BilinearForm::TrialElement : public dolfin::FiniteElement
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
    return 10;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    error("Element is scalar.");
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
    static unsigned int edge_reordering_0[2][2] = {{0, 1}, {1, 0}};
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    int alignment = cell.alignment(1, 0);
    int offset = mesh.topology().size(0);
    nodes[3] = offset + 2*cell.entities(1)[0] + edge_reordering_0[alignment][0];
    nodes[4] = offset + 2*cell.entities(1)[0] + edge_reordering_0[alignment][1];
    alignment = cell.alignment(1, 1);
    nodes[5] = offset + 2*cell.entities(1)[1] + edge_reordering_0[alignment][0];
    nodes[6] = offset + 2*cell.entities(1)[1] + edge_reordering_0[alignment][1];
    alignment = cell.alignment(1, 2);
    nodes[7] = offset + 2*cell.entities(1)[2] + edge_reordering_0[alignment][0];
    nodes[8] = offset + 2*cell.entities(1)[2] + edge_reordering_0[alignment][1];
    offset = offset + 2*mesh.topology().size(1);
    nodes[9] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(6.666666666666667e-01, 3.333333333333333e-01);
    points[4] = map(3.333333333333334e-01, 6.666666666666666e-01);
    points[5] = map(0.000000000000000e+00, 6.666666666666667e-01);
    points[6] = map(0.000000000000000e+00, 3.333333333333334e-01);
    points[7] = map(3.333333333333333e-01, 0.000000000000000e+00);
    points[8] = map(6.666666666666666e-01, 0.000000000000000e+00);
    points[9] = map(3.333333333333333e-01, 3.333333333333333e-01);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 0;
    components[5] = 0;
    components[6] = 0;
    components[7] = 0;
    components[8] = 0;
    components[9] = 0;
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
    FiniteElementSpec s("Lagrange", "triangle", 3);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

BilinearForm::BilinearForm() : dolfin::BilinearForm(0)
{
  // Create finite element for test space
  _test = new TestElement();

  // Create finite element for trial space
  _trial = new TrialElement();
}

// Contribution from the interior
bool BilinearForm::interior_contribution() const { return true; }

void BilinearForm::eval(real block[], const AffineMap& map, real det) const
{
  // Compute geometry tensors
  const real G0_0_0 = det*(map.g00*map.g00 + map.g01*map.g01);
  const real G0_0_1 = det*(map.g00*map.g10 + map.g01*map.g11);
  const real G0_1_0 = det*(map.g10*map.g00 + map.g11*map.g01);
  const real G0_1_1 = det*(map.g10*map.g10 + map.g11*map.g11);

  // Compute element tensor
  block[0] = 4.249999999999996e-01*G0_0_0 + 4.249999999999995e-01*G0_0_1 + 4.249999999999996e-01*G0_1_0 + 4.249999999999996e-01*G0_1_1;
  block[1] = -8.750000000000001e-02*G0_0_0 - 8.750000000000012e-02*G0_1_0;
  block[2] = -8.749999999999999e-02*G0_0_1 - 8.750000000000005e-02*G0_1_1;
  block[3] = -3.749999999999983e-02*G0_0_0 - 3.749999999999998e-02*G0_0_1 - 3.749999999999993e-02*G0_1_0 - 3.750000000000019e-02*G0_1_1;
  block[4] = -3.749999999999992e-02*G0_0_0 - 3.749999999999988e-02*G0_0_1 - 3.749999999999999e-02*G0_1_0 - 3.749999999999987e-02*G0_1_1;
  block[5] = 3.749999999999992e-02*G0_0_0 + 3.374999999999997e-01*G0_0_1 + 3.749999999999998e-02*G0_1_0 + 3.374999999999999e-01*G0_1_1;
  block[6] = 3.750000000000005e-02*G0_0_0 - 6.749999999999992e-01*G0_0_1 + 3.749999999999994e-02*G0_1_0 - 6.749999999999994e-01*G0_1_1;
  block[7] = -6.749999999999996e-01*G0_0_0 + 3.749999999999983e-02*G0_0_1 - 6.749999999999996e-01*G0_1_0 + 3.749999999999981e-02*G0_1_1;
  block[8] = 3.374999999999999e-01*G0_0_0 + 3.750000000000007e-02*G0_0_1 + 3.375000000000000e-01*G0_1_0 + 3.750000000000031e-02*G0_1_1;
  block[9] = 0.000000000000000e+00;
  block[10] = -8.750000000000001e-02*G0_0_0 - 8.750000000000011e-02*G0_0_1;
  block[11] = 4.249999999999998e-01*G0_0_0;
  block[12] = 8.750000000000001e-02*G0_0_1;
  block[13] = 3.750000000000072e-02*G0_0_0 + 7.124999999999999e-01*G0_0_1;
  block[14] = 3.750000000000015e-02*G0_0_0 - 2.999999999999994e-01*G0_0_1;
  block[15] = -3.750000000000015e-02*G0_0_0;
  block[16] = -3.749999999999934e-02*G0_0_0;
  block[17] = 3.374999999999999e-01*G0_0_0 + 3.000000000000002e-01*G0_0_1;
  block[18] = -6.749999999999996e-01*G0_0_0 - 7.124999999999998e-01*G0_0_1;
  block[19] = 0.000000000000000e+00;
  block[20] = -8.750000000000001e-02*G0_1_0 - 8.750000000000005e-02*G0_1_1;
  block[21] = 8.750000000000001e-02*G0_1_0;
  block[22] = 4.249999999999997e-01*G0_1_1;
  block[23] = -2.999999999999994e-01*G0_1_0 + 3.750000000000034e-02*G0_1_1;
  block[24] = 7.124999999999995e-01*G0_1_0 + 3.750000000000017e-02*G0_1_1;
  block[25] = -7.124999999999995e-01*G0_1_0 - 6.749999999999997e-01*G0_1_1;
  block[26] = 2.999999999999999e-01*G0_1_0 + 3.375000000000000e-01*G0_1_1;
  block[27] = -3.749999999999995e-02*G0_1_1;
  block[28] = -3.750000000000028e-02*G0_1_1;
  block[29] = 0.000000000000000e+00;
  block[30] = -3.749999999999984e-02*G0_0_0 - 3.749999999999994e-02*G0_0_1 - 3.749999999999998e-02*G0_1_0 - 3.750000000000021e-02*G0_1_1;
  block[31] = 3.750000000000075e-02*G0_0_0 + 7.124999999999999e-01*G0_1_0;
  block[32] = -2.999999999999994e-01*G0_0_1 + 3.750000000000032e-02*G0_1_1;
  block[33] = 1.687499999999998e+00*G0_0_0 + 8.437500000000000e-01*G0_0_1 + 8.437500000000000e-01*G0_1_0 + 1.687500000000000e+00*G0_1_1;
  block[34] = -3.374999999999996e-01*G0_0_0 + 8.437499999999983e-01*G0_0_1 - 1.687499999999996e-01*G0_1_0 - 3.375000000000002e-01*G0_1_1;
  block[35] = 3.374999999999996e-01*G0_0_0 + 1.687500000000000e-01*G0_0_1 + 1.687499999999996e-01*G0_1_0;
  block[36] = 3.374999999999995e-01*G0_0_0 + 1.687499999999995e-01*G0_0_1 + 1.687500000000005e-01*G0_1_0;
  block[37] = 1.687499999999998e-01*G0_0_1 + 1.687499999999993e-01*G0_1_0 + 3.374999999999996e-01*G0_1_1;
  block[38] = -8.437500000000001e-01*G0_0_1 - 8.437499999999993e-01*G0_1_0 - 1.687500000000000e+00*G0_1_1;
  block[39] = -2.024999999999997e+00*G0_0_0 - 1.012499999999998e+00*G0_0_1 - 1.012500000000001e+00*G0_1_0;
  block[40] = -3.749999999999992e-02*G0_0_0 - 3.749999999999997e-02*G0_0_1 - 3.749999999999989e-02*G0_1_0 - 3.749999999999988e-02*G0_1_1;
  block[41] = 3.750000000000015e-02*G0_0_0 - 2.999999999999994e-01*G0_1_0;
  block[42] = 7.124999999999996e-01*G0_0_1 + 3.750000000000019e-02*G0_1_1;
  block[43] = -3.374999999999996e-01*G0_0_0 - 1.687499999999996e-01*G0_0_1 + 8.437499999999982e-01*G0_1_0 - 3.375000000000002e-01*G0_1_1;
  block[44] = 1.687499999999999e+00*G0_0_0 + 8.437499999999993e-01*G0_0_1 + 8.437499999999991e-01*G0_1_0 + 1.687499999999998e+00*G0_1_1;
  block[45] = -1.687499999999999e+00*G0_0_0 - 8.437499999999997e-01*G0_0_1 - 8.437499999999991e-01*G0_1_0;
  block[46] = 3.375000000000000e-01*G0_0_0 + 1.687500000000003e-01*G0_0_1 + 1.687499999999998e-01*G0_1_0;
  block[47] = 1.687500000000000e-01*G0_0_1 + 1.687499999999997e-01*G0_1_0 + 3.374999999999997e-01*G0_1_1;
  block[48] = 1.687499999999993e-01*G0_0_1 + 1.687499999999997e-01*G0_1_0 + 3.374999999999995e-01*G0_1_1;
  block[49] = -1.012500000000000e+00*G0_0_1 - 1.012499999999998e+00*G0_1_0 - 2.024999999999998e+00*G0_1_1;
  block[50] = 3.749999999999992e-02*G0_0_0 + 3.749999999999997e-02*G0_0_1 + 3.374999999999997e-01*G0_1_0 + 3.374999999999998e-01*G0_1_1;
  block[51] = -3.750000000000013e-02*G0_0_0;
  block[52] = -7.124999999999996e-01*G0_0_1 - 6.749999999999997e-01*G0_1_1;
  block[53] = 3.374999999999996e-01*G0_0_0 + 1.687499999999996e-01*G0_0_1 + 1.687499999999999e-01*G0_1_0;
  block[54] = -1.687499999999999e+00*G0_0_0 - 8.437499999999993e-01*G0_0_1 - 8.437499999999997e-01*G0_1_0;
  block[55] = 1.687499999999999e+00*G0_0_0 + 8.437499999999998e-01*G0_0_1 + 8.437499999999997e-01*G0_1_0 + 1.687499999999999e+00*G0_1_1;
  block[56] = -3.375000000000000e-01*G0_0_0 - 1.687500000000003e-01*G0_0_1 - 1.181249999999999e+00*G0_1_0 - 1.349999999999999e+00*G0_1_1;
  block[57] = -1.687500000000000e-01*G0_0_1 - 1.687500000000001e-01*G0_1_0;
  block[58] = -1.687499999999994e-01*G0_0_1 - 1.687499999999993e-01*G0_1_0;
  block[59] = 1.012500000000000e+00*G0_0_1 + 1.012500000000000e+00*G0_1_0;
  block[60] = 3.750000000000004e-02*G0_0_0 + 3.749999999999994e-02*G0_0_1 - 6.749999999999992e-01*G0_1_0 - 6.749999999999994e-01*G0_1_1;
  block[61] = -3.749999999999933e-02*G0_0_0;
  block[62] = 2.999999999999999e-01*G0_0_1 + 3.375000000000000e-01*G0_1_1;
  block[63] = 3.374999999999995e-01*G0_0_0 + 1.687500000000005e-01*G0_0_1 + 1.687499999999995e-01*G0_1_0;
  block[64] = 3.374999999999999e-01*G0_0_0 + 1.687499999999998e-01*G0_0_1 + 1.687500000000003e-01*G0_1_0;
  block[65] = -3.375000000000000e-01*G0_0_0 - 1.181249999999999e+00*G0_0_1 - 1.687500000000003e-01*G0_1_0 - 1.349999999999999e+00*G0_1_1;
  block[66] = 1.687499999999999e+00*G0_0_0 + 8.437499999999993e-01*G0_0_1 + 8.437499999999993e-01*G0_1_0 + 1.687499999999998e+00*G0_1_1;
  block[67] = 8.437499999999997e-01*G0_0_1 + 8.437499999999994e-01*G0_1_0;
  block[68] = -1.687500000000007e-01*G0_0_1 - 1.687500000000004e-01*G0_1_0;
  block[69] = -2.024999999999999e+00*G0_0_0 - 1.012499999999999e+00*G0_0_1 - 1.012499999999999e+00*G0_1_0;
  block[70] = -6.749999999999996e-01*G0_0_0 - 6.749999999999995e-01*G0_0_1 + 3.749999999999989e-02*G0_1_0 + 3.749999999999981e-02*G0_1_1;
  block[71] = 3.374999999999999e-01*G0_0_0 + 3.000000000000002e-01*G0_1_0;
  block[72] = -3.749999999999995e-02*G0_1_1;
  block[73] = 1.687499999999994e-01*G0_0_1 + 1.687499999999997e-01*G0_1_0 + 3.374999999999996e-01*G0_1_1;
  block[74] = 1.687499999999997e-01*G0_0_1 + 1.687500000000000e-01*G0_1_0 + 3.374999999999997e-01*G0_1_1;
  block[75] = -1.687500000000001e-01*G0_0_1 - 1.687500000000000e-01*G0_1_0;
  block[76] = 8.437499999999996e-01*G0_0_1 + 8.437499999999997e-01*G0_1_0;
  block[77] = 1.687499999999999e+00*G0_0_0 + 8.437500000000001e-01*G0_0_1 + 8.437500000000001e-01*G0_1_0 + 1.687499999999999e+00*G0_1_1;
  block[78] = -1.350000000000000e+00*G0_0_0 - 1.687499999999999e-01*G0_0_1 - 1.181250000000000e+00*G0_1_0 - 3.375000000000002e-01*G0_1_1;
  block[79] = -1.012499999999999e+00*G0_0_1 - 1.012500000000000e+00*G0_1_0 - 2.024999999999999e+00*G0_1_1;
  block[80] = 3.374999999999999e-01*G0_0_0 + 3.374999999999999e-01*G0_0_1 + 3.750000000000007e-02*G0_1_0 + 3.750000000000032e-02*G0_1_1;
  block[81] = -6.749999999999997e-01*G0_0_0 - 7.124999999999998e-01*G0_1_0;
  block[82] = -3.750000000000028e-02*G0_1_1;
  block[83] = -8.437499999999994e-01*G0_0_1 - 8.437500000000001e-01*G0_1_0 - 1.687500000000000e+00*G0_1_1;
  block[84] = 1.687499999999997e-01*G0_0_1 + 1.687499999999993e-01*G0_1_0 + 3.374999999999995e-01*G0_1_1;
  block[85] = -1.687499999999992e-01*G0_0_1 - 1.687499999999994e-01*G0_1_0;
  block[86] = -1.687500000000004e-01*G0_0_1 - 1.687500000000006e-01*G0_1_0;
  block[87] = -1.350000000000000e+00*G0_0_0 - 1.181250000000000e+00*G0_0_1 - 1.687499999999999e-01*G0_1_0 - 3.375000000000002e-01*G0_1_1;
  block[88] = 1.687499999999999e+00*G0_0_0 + 8.437499999999998e-01*G0_0_1 + 8.437499999999998e-01*G0_1_0 + 1.687499999999999e+00*G0_1_1;
  block[89] = 1.012500000000000e+00*G0_0_1 + 1.012500000000001e+00*G0_1_0;
  block[90] = 0.000000000000000e+00;
  block[91] = 0.000000000000000e+00;
  block[92] = 0.000000000000000e+00;
  block[93] = -2.024999999999998e+00*G0_0_0 - 1.012500000000001e+00*G0_0_1 - 1.012499999999998e+00*G0_1_0;
  block[94] = -1.012499999999998e+00*G0_0_1 - 1.012500000000000e+00*G0_1_0 - 2.024999999999998e+00*G0_1_1;
  block[95] = 1.012500000000000e+00*G0_0_1 + 1.012500000000000e+00*G0_1_0;
  block[96] = -2.024999999999999e+00*G0_0_0 - 1.012499999999999e+00*G0_0_1 - 1.012499999999999e+00*G0_1_0;
  block[97] = -1.012500000000000e+00*G0_0_1 - 1.012499999999999e+00*G0_1_0 - 2.024999999999999e+00*G0_1_1;
  block[98] = 1.012500000000001e+00*G0_0_1 + 1.012500000000000e+00*G0_1_0;
  block[99] = 4.049999999999996e+00*G0_0_0 + 2.024999999999998e+00*G0_0_1 + 2.024999999999998e+00*G0_1_0 + 4.049999999999997e+00*G0_1_1;
}

// No contribution from the boundary
bool BilinearForm::boundary_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map, real det, unsigned int facet) const {}

// No contribution from interior boundaries
bool BilinearForm::interior_boundary_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const {}

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class LinearForm : public dolfin::LinearForm
{
public:

  class TestElement;

  class FunctionElement_0;

  LinearForm(Function& w0);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map, real det) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, real det, unsigned int facet) const;

  bool interior_boundary_contribution() const;

  void eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const;

};

class LinearForm::TestElement : public dolfin::FiniteElement
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
    return 10;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    error("Element is scalar.");
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
    static unsigned int edge_reordering_0[2][2] = {{0, 1}, {1, 0}};
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    int alignment = cell.alignment(1, 0);
    int offset = mesh.topology().size(0);
    nodes[3] = offset + 2*cell.entities(1)[0] + edge_reordering_0[alignment][0];
    nodes[4] = offset + 2*cell.entities(1)[0] + edge_reordering_0[alignment][1];
    alignment = cell.alignment(1, 1);
    nodes[5] = offset + 2*cell.entities(1)[1] + edge_reordering_0[alignment][0];
    nodes[6] = offset + 2*cell.entities(1)[1] + edge_reordering_0[alignment][1];
    alignment = cell.alignment(1, 2);
    nodes[7] = offset + 2*cell.entities(1)[2] + edge_reordering_0[alignment][0];
    nodes[8] = offset + 2*cell.entities(1)[2] + edge_reordering_0[alignment][1];
    offset = offset + 2*mesh.topology().size(1);
    nodes[9] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(6.666666666666667e-01, 3.333333333333333e-01);
    points[4] = map(3.333333333333334e-01, 6.666666666666666e-01);
    points[5] = map(0.000000000000000e+00, 6.666666666666667e-01);
    points[6] = map(0.000000000000000e+00, 3.333333333333334e-01);
    points[7] = map(3.333333333333333e-01, 0.000000000000000e+00);
    points[8] = map(6.666666666666666e-01, 0.000000000000000e+00);
    points[9] = map(3.333333333333333e-01, 3.333333333333333e-01);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 0;
    components[5] = 0;
    components[6] = 0;
    components[7] = 0;
    components[8] = 0;
    components[9] = 0;
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
    FiniteElementSpec s("Lagrange", "triangle", 3);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

class LinearForm::FunctionElement_0 : public dolfin::FiniteElement
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
    return 10;
  }

  inline unsigned int shapedim() const
  {
    return 2;
  }

  inline unsigned int tensordim(unsigned int i) const
  {
    error("Element is scalar.");
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
    static unsigned int edge_reordering_0[2][2] = {{0, 1}, {1, 0}};
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    int alignment = cell.alignment(1, 0);
    int offset = mesh.topology().size(0);
    nodes[3] = offset + 2*cell.entities(1)[0] + edge_reordering_0[alignment][0];
    nodes[4] = offset + 2*cell.entities(1)[0] + edge_reordering_0[alignment][1];
    alignment = cell.alignment(1, 1);
    nodes[5] = offset + 2*cell.entities(1)[1] + edge_reordering_0[alignment][0];
    nodes[6] = offset + 2*cell.entities(1)[1] + edge_reordering_0[alignment][1];
    alignment = cell.alignment(1, 2);
    nodes[7] = offset + 2*cell.entities(1)[2] + edge_reordering_0[alignment][0];
    nodes[8] = offset + 2*cell.entities(1)[2] + edge_reordering_0[alignment][1];
    offset = offset + 2*mesh.topology().size(1);
    nodes[9] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(6.666666666666667e-01, 3.333333333333333e-01);
    points[4] = map(3.333333333333334e-01, 6.666666666666666e-01);
    points[5] = map(0.000000000000000e+00, 6.666666666666667e-01);
    points[6] = map(0.000000000000000e+00, 3.333333333333334e-01);
    points[7] = map(3.333333333333333e-01, 0.000000000000000e+00);
    points[8] = map(6.666666666666666e-01, 0.000000000000000e+00);
    points[9] = map(3.333333333333333e-01, 3.333333333333333e-01);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 0;
    components[4] = 0;
    components[5] = 0;
    components[6] = 0;
    components[7] = 0;
    components[8] = 0;
    components[9] = 0;
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
    FiniteElementSpec s("Lagrange", "triangle", 3);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

LinearForm::LinearForm(Function& w0) : dolfin::LinearForm(1)
{
  // Create finite element for test space
  _test = new TestElement();

  // Add functions
  initFunction(0, w0, new FunctionElement_0());
}

// Contribution from the interior
bool LinearForm::interior_contribution() const { return true; }

void LinearForm::eval(real block[], const AffineMap& map, real det) const
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

  // Compute geometry tensors
  const real G0_0 = det*c0_0;
  const real G0_1 = det*c0_1;
  const real G0_2 = det*c0_2;
  const real G0_3 = det*c0_3;
  const real G0_4 = det*c0_4;
  const real G0_5 = det*c0_5;
  const real G0_6 = det*c0_6;
  const real G0_7 = det*c0_7;
  const real G0_8 = det*c0_8;
  const real G0_9 = det*c0_9;

  // Compute element tensor
  block[0] = 5.654761904761898e-03*G0_0 + 8.184523809523801e-04*G0_1 + 8.184523809523793e-04*G0_2 + 2.008928571428568e-03*G0_3 + 2.008928571428570e-03*G0_4 + 1.339285714285713e-03*G0_6 + 1.339285714285709e-03*G0_7 + 2.678571428571425e-03*G0_9;
  block[1] = 8.184523809523800e-04*G0_0 + 5.654761904761894e-03*G0_1 + 8.184523809523791e-04*G0_2 + 1.339285714285704e-03*G0_3 + 2.008928571428569e-03*G0_5 + 2.008928571428569e-03*G0_6 + 1.339285714285709e-03*G0_8 + 2.678571428571413e-03*G0_9;
  block[2] = 8.184523809523793e-04*G0_0 + 8.184523809523792e-04*G0_1 + 5.654761904761899e-03*G0_2 + 1.339285714285713e-03*G0_4 + 1.339285714285709e-03*G0_5 + 2.008928571428571e-03*G0_7 + 2.008928571428570e-03*G0_8 + 2.678571428571422e-03*G0_9;
  block[3] = 2.008928571428568e-03*G0_0 + 1.339285714285704e-03*G0_1 + 4.017857142857136e-02*G0_3 - 1.406249999999998e-02*G0_4 - 1.004464285714285e-02*G0_5 - 4.017857142857142e-03*G0_6 - 1.004464285714285e-02*G0_7 + 2.008928571428571e-02*G0_8 + 1.205357142857144e-02*G0_9;
  block[4] = 2.008928571428570e-03*G0_0 + 1.339285714285713e-03*G0_2 - 1.406249999999998e-02*G0_3 + 4.017857142857139e-02*G0_4 + 2.008928571428570e-02*G0_5 - 1.004464285714284e-02*G0_6 - 4.017857142857136e-03*G0_7 - 1.004464285714285e-02*G0_8 + 1.205357142857142e-02*G0_9;
  block[5] = 2.008928571428569e-03*G0_1 + 1.339285714285710e-03*G0_2 - 1.004464285714285e-02*G0_3 + 2.008928571428570e-02*G0_4 + 4.017857142857139e-02*G0_5 - 1.406249999999999e-02*G0_6 - 1.004464285714284e-02*G0_7 - 4.017857142857136e-03*G0_8 + 1.205357142857142e-02*G0_9;
  block[6] = 1.339285714285712e-03*G0_0 + 2.008928571428569e-03*G0_1 - 4.017857142857141e-03*G0_3 - 1.004464285714284e-02*G0_4 - 1.406249999999999e-02*G0_5 + 4.017857142857140e-02*G0_6 + 2.008928571428570e-02*G0_7 - 1.004464285714284e-02*G0_8 + 1.205357142857144e-02*G0_9;
  block[7] = 1.339285714285709e-03*G0_0 + 2.008928571428571e-03*G0_2 - 1.004464285714285e-02*G0_3 - 4.017857142857136e-03*G0_4 - 1.004464285714284e-02*G0_5 + 2.008928571428569e-02*G0_6 + 4.017857142857141e-02*G0_7 - 1.406249999999998e-02*G0_8 + 1.205357142857142e-02*G0_9;
  block[8] = 1.339285714285710e-03*G0_1 + 2.008928571428570e-03*G0_2 + 2.008928571428570e-02*G0_3 - 1.004464285714284e-02*G0_4 - 4.017857142857136e-03*G0_5 - 1.004464285714284e-02*G0_6 - 1.406249999999999e-02*G0_7 + 4.017857142857141e-02*G0_8 + 1.205357142857144e-02*G0_9;
  block[9] = 2.678571428571425e-03*G0_0 + 2.678571428571414e-03*G0_1 + 2.678571428571422e-03*G0_2 + 1.205357142857144e-02*G0_3 + 1.205357142857142e-02*G0_4 + 1.205357142857142e-02*G0_5 + 1.205357142857144e-02*G0_6 + 1.205357142857142e-02*G0_7 + 1.205357142857145e-02*G0_8 + 1.446428571428571e-01*G0_9;
}

// No contribution from the boundary
bool LinearForm::boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map, real det, unsigned int facet) const {}

// No contribution from interior boundaries
bool LinearForm::interior_boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const {}

} }

#endif
