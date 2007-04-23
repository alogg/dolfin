// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.5.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU LGPL Version 2.1.

#ifndef __POISSON3D_2_H
#define __POISSON3D_2_H

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

namespace dolfin { namespace Poisson3D_2 {

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
    return 3;
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
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    nodes[3] = cell.entities(0)[3];
    int offset = mesh.topology().size(0);
    nodes[4] = offset + cell.entities(1)[0];
    nodes[5] = offset + cell.entities(1)[1];
    nodes[6] = offset + cell.entities(1)[2];
    nodes[7] = offset + cell.entities(1)[3];
    nodes[8] = offset + cell.entities(1)[4];
    nodes[9] = offset + cell.entities(1)[5];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[4] = map(5.000000000000000e-01, 5.000000000000000e-01, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 5.000000000000000e-01, 0.000000000000000e+00);
    points[6] = map(5.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
    points[7] = map(0.000000000000000e+00, 0.000000000000000e+00, 5.000000000000000e-01);
    points[8] = map(5.000000000000000e-01, 0.000000000000000e+00, 5.000000000000000e-01);
    points[9] = map(0.000000000000000e+00, 5.000000000000000e-01, 5.000000000000000e-01);
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
    FiniteElementSpec s("Lagrange", "tetrahedron", 2);
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
    return 3;
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
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    nodes[3] = cell.entities(0)[3];
    int offset = mesh.topology().size(0);
    nodes[4] = offset + cell.entities(1)[0];
    nodes[5] = offset + cell.entities(1)[1];
    nodes[6] = offset + cell.entities(1)[2];
    nodes[7] = offset + cell.entities(1)[3];
    nodes[8] = offset + cell.entities(1)[4];
    nodes[9] = offset + cell.entities(1)[5];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[4] = map(5.000000000000000e-01, 5.000000000000000e-01, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 5.000000000000000e-01, 0.000000000000000e+00);
    points[6] = map(5.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
    points[7] = map(0.000000000000000e+00, 0.000000000000000e+00, 5.000000000000000e-01);
    points[8] = map(5.000000000000000e-01, 0.000000000000000e+00, 5.000000000000000e-01);
    points[9] = map(0.000000000000000e+00, 5.000000000000000e-01, 5.000000000000000e-01);
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
    FiniteElementSpec s("Lagrange", "tetrahedron", 2);
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
  const real G0_0_0 = det*(map.g00*map.g00 + map.g01*map.g01 + map.g02*map.g02);
  const real G0_0_1 = det*(map.g00*map.g10 + map.g01*map.g11 + map.g02*map.g12);
  const real G0_0_2 = det*(map.g00*map.g20 + map.g01*map.g21 + map.g02*map.g22);
  const real G0_1_0 = det*(map.g10*map.g00 + map.g11*map.g01 + map.g12*map.g02);
  const real G0_1_1 = det*(map.g10*map.g10 + map.g11*map.g11 + map.g12*map.g12);
  const real G0_1_2 = det*(map.g10*map.g20 + map.g11*map.g21 + map.g12*map.g22);
  const real G0_2_0 = det*(map.g20*map.g00 + map.g21*map.g01 + map.g22*map.g02);
  const real G0_2_1 = det*(map.g20*map.g10 + map.g21*map.g11 + map.g22*map.g12);
  const real G0_2_2 = det*(map.g20*map.g20 + map.g21*map.g21 + map.g22*map.g22);

  // Compute element tensor
  block[0] = 9.999999999999978e-02*G0_0_0 + 9.999999999999977e-02*G0_0_1 + 9.999999999999977e-02*G0_0_2 + 9.999999999999977e-02*G0_1_0 + 9.999999999999974e-02*G0_1_1 + 9.999999999999976e-02*G0_1_2 + 9.999999999999977e-02*G0_2_0 + 9.999999999999974e-02*G0_2_1 + 9.999999999999976e-02*G0_2_2;
  block[1] = 3.333333333333328e-02*G0_0_0 + 3.333333333333326e-02*G0_1_0 + 3.333333333333327e-02*G0_2_0;
  block[2] = 3.333333333333325e-02*G0_0_1 + 3.333333333333325e-02*G0_1_1 + 3.333333333333324e-02*G0_2_1;
  block[3] = 3.333333333333323e-02*G0_0_2 + 3.333333333333324e-02*G0_1_2 + 3.333333333333323e-02*G0_2_2;
  block[4] = 3.333333333333326e-02*G0_0_0 + 3.333333333333328e-02*G0_0_1 + 3.333333333333324e-02*G0_1_0 + 3.333333333333326e-02*G0_1_1 + 3.333333333333324e-02*G0_2_0 + 3.333333333333328e-02*G0_2_1;
  block[5] = -3.333333333333328e-02*G0_0_0 - 1.333333333333330e-01*G0_0_1 - 3.333333333333328e-02*G0_0_2 - 3.333333333333327e-02*G0_1_0 - 1.333333333333329e-01*G0_1_1 - 3.333333333333326e-02*G0_1_2 - 3.333333333333326e-02*G0_2_0 - 1.333333333333329e-01*G0_2_1 - 3.333333333333326e-02*G0_2_2;
  block[6] = -1.333333333333331e-01*G0_0_0 - 3.333333333333328e-02*G0_0_1 - 3.333333333333330e-02*G0_0_2 - 1.333333333333330e-01*G0_1_0 - 3.333333333333325e-02*G0_1_1 - 3.333333333333326e-02*G0_1_2 - 1.333333333333331e-01*G0_2_0 - 3.333333333333327e-02*G0_2_1 - 3.333333333333328e-02*G0_2_2;
  block[7] = -3.333333333333326e-02*G0_0_0 - 3.333333333333326e-02*G0_0_1 - 1.333333333333330e-01*G0_0_2 - 3.333333333333324e-02*G0_1_0 - 3.333333333333324e-02*G0_1_1 - 1.333333333333330e-01*G0_1_2 - 3.333333333333324e-02*G0_2_0 - 3.333333333333324e-02*G0_2_1 - 1.333333333333330e-01*G0_2_2;
  block[8] = 3.333333333333325e-02*G0_0_0 + 3.333333333333329e-02*G0_0_2 + 3.333333333333324e-02*G0_1_0 + 3.333333333333326e-02*G0_1_2 + 3.333333333333324e-02*G0_2_0 + 3.333333333333328e-02*G0_2_2;
  block[9] = 3.333333333333325e-02*G0_0_1 + 3.333333333333326e-02*G0_0_2 + 3.333333333333324e-02*G0_1_1 + 3.333333333333326e-02*G0_1_2 + 3.333333333333324e-02*G0_2_1 + 3.333333333333326e-02*G0_2_2;
  block[10] = 3.333333333333328e-02*G0_0_0 + 3.333333333333326e-02*G0_0_1 + 3.333333333333327e-02*G0_0_2;
  block[11] = 9.999999999999976e-02*G0_0_0;
  block[12] = -3.333333333333321e-02*G0_0_1;
  block[13] = -3.333333333333324e-02*G0_0_2;
  block[14] = -3.333333333333324e-02*G0_0_0 + 9.999999999999974e-02*G0_0_1;
  block[15] = 3.333333333333324e-02*G0_0_0 + 3.333333333333324e-02*G0_0_2;
  block[16] = -1.333333333333331e-01*G0_0_0 - 9.999999999999973e-02*G0_0_1 - 9.999999999999974e-02*G0_0_2;
  block[17] = 3.333333333333328e-02*G0_0_0 + 3.333333333333326e-02*G0_0_1;
  block[18] = -3.333333333333328e-02*G0_0_0 + 9.999999999999976e-02*G0_0_2;
  block[19] = -3.333333333333327e-02*G0_0_1 - 3.333333333333324e-02*G0_0_2;
  block[20] = 3.333333333333324e-02*G0_1_0 + 3.333333333333325e-02*G0_1_1 + 3.333333333333324e-02*G0_1_2;
  block[21] = -3.333333333333321e-02*G0_1_0;
  block[22] = 9.999999999999969e-02*G0_1_1;
  block[23] = -3.333333333333323e-02*G0_1_2;
  block[24] = 9.999999999999970e-02*G0_1_0 - 3.333333333333321e-02*G0_1_1;
  block[25] = -9.999999999999971e-02*G0_1_0 - 1.333333333333329e-01*G0_1_1 - 9.999999999999970e-02*G0_1_2;
  block[26] = 3.333333333333321e-02*G0_1_1 + 3.333333333333321e-02*G0_1_2;
  block[27] = 3.333333333333321e-02*G0_1_0 + 3.333333333333321e-02*G0_1_1;
  block[28] = -3.333333333333321e-02*G0_1_0 - 3.333333333333322e-02*G0_1_2;
  block[29] = -3.333333333333321e-02*G0_1_1 + 9.999999999999971e-02*G0_1_2;
  block[30] = 3.333333333333323e-02*G0_2_0 + 3.333333333333324e-02*G0_2_1 + 3.333333333333322e-02*G0_2_2;
  block[31] = -3.333333333333324e-02*G0_2_0;
  block[32] = -3.333333333333323e-02*G0_2_1;
  block[33] = 9.999999999999970e-02*G0_2_2;
  block[34] = -3.333333333333317e-02*G0_2_0 - 3.333333333333320e-02*G0_2_1;
  block[35] = 3.333333333333319e-02*G0_2_0 + 3.333333333333317e-02*G0_2_2;
  block[36] = 3.333333333333319e-02*G0_2_1 + 3.333333333333318e-02*G0_2_2;
  block[37] = -9.999999999999977e-02*G0_2_0 - 9.999999999999977e-02*G0_2_1 - 1.333333333333329e-01*G0_2_2;
  block[38] = 9.999999999999977e-02*G0_2_0 - 3.333333333333320e-02*G0_2_2;
  block[39] = 9.999999999999976e-02*G0_2_1 - 3.333333333333317e-02*G0_2_2;
  block[40] = 3.333333333333326e-02*G0_0_0 + 3.333333333333324e-02*G0_0_1 + 3.333333333333324e-02*G0_0_2 + 3.333333333333328e-02*G0_1_0 + 3.333333333333326e-02*G0_1_1 + 3.333333333333328e-02*G0_1_2;
  block[41] = -3.333333333333324e-02*G0_0_0 + 9.999999999999974e-02*G0_1_0;
  block[42] = 9.999999999999970e-02*G0_0_1 - 3.333333333333321e-02*G0_1_1;
  block[43] = -3.333333333333317e-02*G0_0_2 - 3.333333333333320e-02*G0_1_2;
  block[44] = 2.666666666666659e-01*G0_0_0 + 1.333333333333329e-01*G0_0_1 + 1.333333333333329e-01*G0_1_0 + 2.666666666666659e-01*G0_1_1;
  block[45] = -2.666666666666660e-01*G0_0_0 - 1.333333333333330e-01*G0_0_1 - 2.666666666666659e-01*G0_0_2 - 1.333333333333330e-01*G0_1_0 - 1.333333333333330e-01*G0_1_2;
  block[46] = -1.333333333333329e-01*G0_0_1 - 1.333333333333329e-01*G0_0_2 - 1.333333333333330e-01*G0_1_0 - 2.666666666666658e-01*G0_1_1 - 2.666666666666658e-01*G0_1_2;
  block[47] = -1.333333333333331e-01*G0_0_0 - 1.333333333333330e-01*G0_0_1 - 1.333333333333330e-01*G0_1_0 - 1.333333333333330e-01*G0_1_1;
  block[48] = 1.333333333333331e-01*G0_0_0 + 1.333333333333330e-01*G0_0_2 + 1.333333333333330e-01*G0_1_0 + 2.666666666666659e-01*G0_1_2;
  block[49] = 1.333333333333330e-01*G0_0_1 + 2.666666666666659e-01*G0_0_2 + 1.333333333333330e-01*G0_1_1 + 1.333333333333329e-01*G0_1_2;
  block[50] = -3.333333333333328e-02*G0_0_0 - 3.333333333333327e-02*G0_0_1 - 3.333333333333326e-02*G0_0_2 - 1.333333333333330e-01*G0_1_0 - 1.333333333333330e-01*G0_1_1 - 1.333333333333330e-01*G0_1_2 - 3.333333333333328e-02*G0_2_0 - 3.333333333333326e-02*G0_2_1 - 3.333333333333326e-02*G0_2_2;
  block[51] = 3.333333333333325e-02*G0_0_0 + 3.333333333333324e-02*G0_2_0;
  block[52] = -9.999999999999973e-02*G0_0_1 - 1.333333333333329e-01*G0_1_1 - 9.999999999999970e-02*G0_2_1;
  block[53] = 3.333333333333319e-02*G0_0_2 + 3.333333333333317e-02*G0_2_2;
  block[54] = -2.666666666666660e-01*G0_0_0 - 1.333333333333330e-01*G0_0_1 - 1.333333333333330e-01*G0_1_0 - 2.666666666666659e-01*G0_2_0 - 1.333333333333329e-01*G0_2_1;
  block[55] = 2.666666666666661e-01*G0_0_0 + 1.333333333333330e-01*G0_0_1 + 2.666666666666660e-01*G0_0_2 + 1.333333333333330e-01*G0_1_0 + 2.666666666666658e-01*G0_1_1 + 1.333333333333330e-01*G0_1_2 + 2.666666666666660e-01*G0_2_0 + 1.333333333333330e-01*G0_2_1 + 2.666666666666659e-01*G0_2_2;
  block[56] = 1.333333333333330e-01*G0_0_1 + 1.333333333333329e-01*G0_0_2 + 1.333333333333331e-01*G0_1_0 + 1.333333333333329e-01*G0_2_1 + 1.333333333333329e-01*G0_2_2;
  block[57] = 1.333333333333331e-01*G0_0_0 + 1.333333333333331e-01*G0_0_1 + 1.333333333333330e-01*G0_1_2 + 1.333333333333331e-01*G0_2_0 + 1.333333333333331e-01*G0_2_1;
  block[58] = -1.333333333333331e-01*G0_0_0 - 1.333333333333330e-01*G0_0_2 - 1.333333333333331e-01*G0_2_0 - 1.333333333333330e-01*G0_2_2;
  block[59] = -1.333333333333331e-01*G0_0_1 - 2.666666666666660e-01*G0_0_2 - 1.333333333333330e-01*G0_1_2 - 1.333333333333331e-01*G0_2_1 - 2.666666666666659e-01*G0_2_2;
  block[60] = -1.333333333333331e-01*G0_0_0 - 1.333333333333330e-01*G0_0_1 - 1.333333333333330e-01*G0_0_2 - 3.333333333333330e-02*G0_1_0 - 3.333333333333325e-02*G0_1_1 - 3.333333333333327e-02*G0_1_2 - 3.333333333333330e-02*G0_2_0 - 3.333333333333326e-02*G0_2_1 - 3.333333333333328e-02*G0_2_2;
  block[61] = -1.333333333333331e-01*G0_0_0 - 9.999999999999973e-02*G0_1_0 - 9.999999999999974e-02*G0_2_0;
  block[62] = 3.333333333333321e-02*G0_1_1 + 3.333333333333322e-02*G0_2_1;
  block[63] = 3.333333333333319e-02*G0_1_2 + 3.333333333333319e-02*G0_2_2;
  block[64] = -1.333333333333330e-01*G0_0_1 - 1.333333333333329e-01*G0_1_0 - 2.666666666666659e-01*G0_1_1 - 1.333333333333329e-01*G0_2_0 - 2.666666666666659e-01*G0_2_1;
  block[65] = 1.333333333333330e-01*G0_0_1 + 1.333333333333330e-01*G0_1_0 + 1.333333333333329e-01*G0_1_2 + 1.333333333333329e-01*G0_2_0 + 1.333333333333329e-01*G0_2_2;
  block[66] = 2.666666666666662e-01*G0_0_0 + 1.333333333333330e-01*G0_0_1 + 1.333333333333330e-01*G0_0_2 + 1.333333333333330e-01*G0_1_0 + 2.666666666666658e-01*G0_1_1 + 2.666666666666658e-01*G0_1_2 + 1.333333333333330e-01*G0_2_0 + 2.666666666666658e-01*G0_2_1 + 2.666666666666658e-01*G0_2_2;
  block[67] = 1.333333333333330e-01*G0_0_2 + 1.333333333333330e-01*G0_1_0 + 1.333333333333329e-01*G0_1_1 + 1.333333333333330e-01*G0_2_0 + 1.333333333333329e-01*G0_2_1;
  block[68] = -1.333333333333331e-01*G0_0_2 - 1.333333333333330e-01*G0_1_0 - 2.666666666666659e-01*G0_1_2 - 1.333333333333330e-01*G0_2_0 - 2.666666666666659e-01*G0_2_2;
  block[69] = -1.333333333333329e-01*G0_1_1 - 1.333333333333329e-01*G0_1_2 - 1.333333333333329e-01*G0_2_1 - 1.333333333333329e-01*G0_2_2;
  block[70] = -3.333333333333326e-02*G0_0_0 - 3.333333333333324e-02*G0_0_1 - 3.333333333333324e-02*G0_0_2 - 3.333333333333325e-02*G0_1_0 - 3.333333333333324e-02*G0_1_1 - 3.333333333333324e-02*G0_1_2 - 1.333333333333330e-01*G0_2_0 - 1.333333333333329e-01*G0_2_1 - 1.333333333333329e-01*G0_2_2;
  block[71] = 3.333333333333328e-02*G0_0_0 + 3.333333333333327e-02*G0_1_0;
  block[72] = 3.333333333333321e-02*G0_0_1 + 3.333333333333321e-02*G0_1_1;
  block[73] = -9.999999999999977e-02*G0_0_2 - 9.999999999999977e-02*G0_1_2 - 1.333333333333329e-01*G0_2_2;
  block[74] = -1.333333333333331e-01*G0_0_0 - 1.333333333333330e-01*G0_0_1 - 1.333333333333330e-01*G0_1_0 - 1.333333333333330e-01*G0_1_1;
  block[75] = 1.333333333333331e-01*G0_0_0 + 1.333333333333331e-01*G0_0_2 + 1.333333333333331e-01*G0_1_0 + 1.333333333333331e-01*G0_1_2 + 1.333333333333329e-01*G0_2_1;
  block[76] = 1.333333333333330e-01*G0_0_1 + 1.333333333333330e-01*G0_0_2 + 1.333333333333329e-01*G0_1_1 + 1.333333333333329e-01*G0_1_2 + 1.333333333333330e-01*G0_2_0;
  block[77] = 2.666666666666662e-01*G0_0_0 + 2.666666666666661e-01*G0_0_1 + 1.333333333333331e-01*G0_0_2 + 2.666666666666661e-01*G0_1_0 + 2.666666666666661e-01*G0_1_1 + 1.333333333333331e-01*G0_1_2 + 1.333333333333331e-01*G0_2_0 + 1.333333333333331e-01*G0_2_1 + 2.666666666666659e-01*G0_2_2;
  block[78] = -2.666666666666662e-01*G0_0_0 - 1.333333333333330e-01*G0_0_2 - 2.666666666666661e-01*G0_1_0 - 1.333333333333330e-01*G0_1_2 - 1.333333333333331e-01*G0_2_0;
  block[79] = -2.666666666666661e-01*G0_0_1 - 1.333333333333331e-01*G0_0_2 - 2.666666666666661e-01*G0_1_1 - 1.333333333333330e-01*G0_1_2 - 1.333333333333330e-01*G0_2_1;
  block[80] = 3.333333333333324e-02*G0_0_0 + 3.333333333333324e-02*G0_0_1 + 3.333333333333323e-02*G0_0_2 + 3.333333333333330e-02*G0_2_0 + 3.333333333333326e-02*G0_2_1 + 3.333333333333328e-02*G0_2_2;
  block[81] = -3.333333333333328e-02*G0_0_0 + 9.999999999999976e-02*G0_2_0;
  block[82] = -3.333333333333321e-02*G0_0_1 - 3.333333333333321e-02*G0_2_1;
  block[83] = 9.999999999999977e-02*G0_0_2 - 3.333333333333319e-02*G0_2_2;
  block[84] = 1.333333333333331e-01*G0_0_0 + 1.333333333333330e-01*G0_0_1 + 1.333333333333330e-01*G0_2_0 + 2.666666666666659e-01*G0_2_1;
  block[85] = -1.333333333333331e-01*G0_0_0 - 1.333333333333331e-01*G0_0_2 - 1.333333333333331e-01*G0_2_0 - 1.333333333333330e-01*G0_2_2;
  block[86] = -1.333333333333330e-01*G0_0_1 - 1.333333333333330e-01*G0_0_2 - 1.333333333333331e-01*G0_2_0 - 2.666666666666659e-01*G0_2_1 - 2.666666666666659e-01*G0_2_2;
  block[87] = -2.666666666666662e-01*G0_0_0 - 2.666666666666661e-01*G0_0_1 - 1.333333333333331e-01*G0_0_2 - 1.333333333333331e-01*G0_2_0 - 1.333333333333330e-01*G0_2_1;
  block[88] = 2.666666666666662e-01*G0_0_0 + 1.333333333333331e-01*G0_0_2 + 1.333333333333331e-01*G0_2_0 + 2.666666666666660e-01*G0_2_2;
  block[89] = 2.666666666666661e-01*G0_0_1 + 1.333333333333331e-01*G0_0_2 + 1.333333333333330e-01*G0_2_1 + 1.333333333333330e-01*G0_2_2;
  block[90] = 3.333333333333324e-02*G0_1_0 + 3.333333333333324e-02*G0_1_1 + 3.333333333333323e-02*G0_1_2 + 3.333333333333326e-02*G0_2_0 + 3.333333333333326e-02*G0_2_1 + 3.333333333333326e-02*G0_2_2;
  block[91] = -3.333333333333328e-02*G0_1_0 - 3.333333333333324e-02*G0_2_0;
  block[92] = -3.333333333333321e-02*G0_1_1 + 9.999999999999971e-02*G0_2_1;
  block[93] = 9.999999999999973e-02*G0_1_2 - 3.333333333333317e-02*G0_2_2;
  block[94] = 1.333333333333330e-01*G0_1_0 + 1.333333333333330e-01*G0_1_1 + 2.666666666666659e-01*G0_2_0 + 1.333333333333329e-01*G0_2_1;
  block[95] = -1.333333333333331e-01*G0_1_0 - 1.333333333333331e-01*G0_1_2 - 2.666666666666659e-01*G0_2_0 - 1.333333333333330e-01*G0_2_1 - 2.666666666666659e-01*G0_2_2;
  block[96] = -1.333333333333329e-01*G0_1_1 - 1.333333333333329e-01*G0_1_2 - 1.333333333333329e-01*G0_2_1 - 1.333333333333329e-01*G0_2_2;
  block[97] = -2.666666666666661e-01*G0_1_0 - 2.666666666666661e-01*G0_1_1 - 1.333333333333330e-01*G0_1_2 - 1.333333333333331e-01*G0_2_0 - 1.333333333333330e-01*G0_2_1;
  block[98] = 2.666666666666661e-01*G0_1_0 + 1.333333333333330e-01*G0_1_2 + 1.333333333333331e-01*G0_2_0 + 1.333333333333330e-01*G0_2_2;
  block[99] = 2.666666666666660e-01*G0_1_1 + 1.333333333333330e-01*G0_1_2 + 1.333333333333330e-01*G0_2_1 + 2.666666666666659e-01*G0_2_2;
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
    return 3;
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
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    nodes[3] = cell.entities(0)[3];
    int offset = mesh.topology().size(0);
    nodes[4] = offset + cell.entities(1)[0];
    nodes[5] = offset + cell.entities(1)[1];
    nodes[6] = offset + cell.entities(1)[2];
    nodes[7] = offset + cell.entities(1)[3];
    nodes[8] = offset + cell.entities(1)[4];
    nodes[9] = offset + cell.entities(1)[5];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[4] = map(5.000000000000000e-01, 5.000000000000000e-01, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 5.000000000000000e-01, 0.000000000000000e+00);
    points[6] = map(5.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
    points[7] = map(0.000000000000000e+00, 0.000000000000000e+00, 5.000000000000000e-01);
    points[8] = map(5.000000000000000e-01, 0.000000000000000e+00, 5.000000000000000e-01);
    points[9] = map(0.000000000000000e+00, 5.000000000000000e-01, 5.000000000000000e-01);
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
    FiniteElementSpec s("Lagrange", "tetrahedron", 2);
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
    return 3;
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
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    nodes[3] = cell.entities(0)[3];
    int offset = mesh.topology().size(0);
    nodes[4] = offset + cell.entities(1)[0];
    nodes[5] = offset + cell.entities(1)[1];
    nodes[6] = offset + cell.entities(1)[2];
    nodes[7] = offset + cell.entities(1)[3];
    nodes[8] = offset + cell.entities(1)[4];
    nodes[9] = offset + cell.entities(1)[5];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[4] = map(5.000000000000000e-01, 5.000000000000000e-01, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 5.000000000000000e-01, 0.000000000000000e+00);
    points[6] = map(5.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
    points[7] = map(0.000000000000000e+00, 0.000000000000000e+00, 5.000000000000000e-01);
    points[8] = map(5.000000000000000e-01, 0.000000000000000e+00, 5.000000000000000e-01);
    points[9] = map(0.000000000000000e+00, 5.000000000000000e-01, 5.000000000000000e-01);
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
    FiniteElementSpec s("Lagrange", "tetrahedron", 2);
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
  block[0] = 2.380952380952376e-03*G0_0 + 3.968253968253965e-04*G0_1 + 3.968253968253968e-04*G0_2 + 3.968253968253969e-04*G0_3 - 2.380952380952378e-03*G0_4 - 1.587301587301586e-03*G0_5 - 1.587301587301588e-03*G0_6 - 1.587301587301587e-03*G0_7 - 2.380952380952379e-03*G0_8 - 2.380952380952379e-03*G0_9;
  block[1] = 3.968253968253966e-04*G0_0 + 2.380952380952378e-03*G0_1 + 3.968253968253959e-04*G0_2 + 3.968253968253961e-04*G0_3 - 1.587301587301583e-03*G0_4 - 2.380952380952376e-03*G0_5 - 1.587301587301585e-03*G0_6 - 2.380952380952378e-03*G0_7 - 1.587301587301585e-03*G0_8 - 2.380952380952378e-03*G0_9;
  block[2] = 3.968253968253968e-04*G0_0 + 3.968253968253958e-04*G0_1 + 2.380952380952378e-03*G0_2 + 3.968253968253968e-04*G0_3 - 1.587301587301586e-03*G0_4 - 1.587301587301585e-03*G0_5 - 2.380952380952379e-03*G0_6 - 2.380952380952379e-03*G0_7 - 2.380952380952381e-03*G0_8 - 1.587301587301585e-03*G0_9;
  block[3] = 3.968253968253969e-04*G0_0 + 3.968253968253961e-04*G0_1 + 3.968253968253968e-04*G0_2 + 2.380952380952378e-03*G0_3 - 2.380952380952379e-03*G0_4 - 2.380952380952379e-03*G0_5 - 2.380952380952379e-03*G0_6 - 1.587301587301586e-03*G0_7 - 1.587301587301587e-03*G0_8 - 1.587301587301587e-03*G0_9;
  block[4] = -2.380952380952378e-03*G0_0 - 1.587301587301582e-03*G0_1 - 1.587301587301585e-03*G0_2 - 2.380952380952379e-03*G0_3 + 1.269841269841269e-02*G0_4 + 6.349206349206338e-03*G0_5 + 6.349206349206346e-03*G0_6 + 3.174603174603171e-03*G0_7 + 6.349206349206349e-03*G0_8 + 6.349206349206347e-03*G0_9;
  block[5] = -1.587301587301586e-03*G0_0 - 2.380952380952376e-03*G0_1 - 1.587301587301585e-03*G0_2 - 2.380952380952379e-03*G0_3 + 6.349206349206338e-03*G0_4 + 1.269841269841268e-02*G0_5 + 6.349206349206347e-03*G0_6 + 6.349206349206343e-03*G0_7 + 3.174603174603173e-03*G0_8 + 6.349206349206341e-03*G0_9;
  block[6] = -1.587301587301587e-03*G0_0 - 1.587301587301585e-03*G0_1 - 2.380952380952379e-03*G0_2 - 2.380952380952379e-03*G0_3 + 6.349206349206346e-03*G0_4 + 6.349206349206347e-03*G0_5 + 1.269841269841268e-02*G0_6 + 6.349206349206343e-03*G0_7 + 6.349206349206345e-03*G0_8 + 3.174603174603175e-03*G0_9;
  block[7] = -1.587301587301587e-03*G0_0 - 2.380952380952378e-03*G0_1 - 2.380952380952379e-03*G0_2 - 1.587301587301586e-03*G0_3 + 3.174603174603172e-03*G0_4 + 6.349206349206343e-03*G0_5 + 6.349206349206343e-03*G0_6 + 1.269841269841268e-02*G0_7 + 6.349206349206342e-03*G0_8 + 6.349206349206343e-03*G0_9;
  block[8] = -2.380952380952379e-03*G0_0 - 1.587301587301584e-03*G0_1 - 2.380952380952381e-03*G0_2 - 1.587301587301587e-03*G0_3 + 6.349206349206349e-03*G0_4 + 3.174603174603173e-03*G0_5 + 6.349206349206345e-03*G0_6 + 6.349206349206342e-03*G0_7 + 1.269841269841270e-02*G0_8 + 6.349206349206348e-03*G0_9;
  block[9] = -2.380952380952379e-03*G0_0 - 2.380952380952377e-03*G0_1 - 1.587301587301586e-03*G0_2 - 1.587301587301587e-03*G0_3 + 6.349206349206346e-03*G0_4 + 6.349206349206341e-03*G0_5 + 3.174603174603175e-03*G0_6 + 6.349206349206343e-03*G0_7 + 6.349206349206347e-03*G0_8 + 1.269841269841269e-02*G0_9;
}

// No contribution from the boundary
bool LinearForm::boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map, real det, unsigned int facet) const {}

// No contribution from interior boundaries
bool LinearForm::interior_boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const {}

} }

#endif
