// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.5.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __VECTORPOISSON_H
#define __VECTORPOISSON_H

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

namespace dolfin { namespace VectorPoisson {

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
    tensordims = new unsigned int [1];
    tensordims[0] = 2;

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
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    int offset = mesh.topology().size(0);
    nodes[3] = offset + cell.entities(0)[0];
    nodes[4] = offset + cell.entities(0)[1];
    nodes[5] = offset + cell.entities(0)[2];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[4] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 1;
    components[4] = 1;
    components[5] = 1;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0);
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
    FiniteElementSpec s("Vector Lagrange", "triangle", 1, 2);
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
    tensordims = new unsigned int [1];
    tensordims[0] = 2;

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
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    int offset = mesh.topology().size(0);
    nodes[3] = offset + cell.entities(0)[0];
    nodes[4] = offset + cell.entities(0)[1];
    nodes[5] = offset + cell.entities(0)[2];
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00);
    points[4] = map(1.000000000000000e+00, 0.000000000000000e+00);
    points[5] = map(0.000000000000000e+00, 1.000000000000000e+00);
    components[0] = 0;
    components[1] = 0;
    components[2] = 0;
    components[3] = 1;
    components[4] = 1;
    components[5] = 1;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0);
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
    FiniteElementSpec s("Vector Lagrange", "triangle", 1, 2);
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
  const real G0_0_0 = det*map.g00*map.g00 + det*map.g00*map.g00 + det*map.g00*map.g00 + det*map.g00*map.g00 + det*map.g01*map.g01 + det*map.g01*map.g01;
  const real G0_0_1 = det*map.g00*map.g10 + det*map.g00*map.g10 + det*map.g00*map.g10 + det*map.g00*map.g10 + det*map.g01*map.g11 + det*map.g01*map.g11;
  const real G0_1_0 = det*map.g10*map.g00 + det*map.g10*map.g00 + det*map.g10*map.g00 + det*map.g10*map.g00 + det*map.g11*map.g01 + det*map.g11*map.g01;
  const real G0_1_1 = det*map.g10*map.g10 + det*map.g10*map.g10 + det*map.g10*map.g10 + det*map.g10*map.g10 + det*map.g11*map.g11 + det*map.g11*map.g11;
  const real G1_0_0 = det*map.g01*map.g00 + det*map.g01*map.g00;
  const real G1_0_1 = det*map.g01*map.g10 + det*map.g01*map.g10;
  const real G1_1_0 = det*map.g11*map.g00 + det*map.g11*map.g00;
  const real G1_1_1 = det*map.g11*map.g10 + det*map.g11*map.g10;
  const real G2_0_0 = det*map.g00*map.g01 + det*map.g00*map.g01;
  const real G2_0_1 = det*map.g00*map.g11 + det*map.g00*map.g11;
  const real G2_1_0 = det*map.g10*map.g01 + det*map.g10*map.g01;
  const real G2_1_1 = det*map.g10*map.g11 + det*map.g10*map.g11;
  const real G3_0_0 = det*map.g00*map.g00 + det*map.g00*map.g00 + det*map.g01*map.g01 + det*map.g01*map.g01 + det*map.g01*map.g01 + det*map.g01*map.g01;
  const real G3_0_1 = det*map.g00*map.g10 + det*map.g00*map.g10 + det*map.g01*map.g11 + det*map.g01*map.g11 + det*map.g01*map.g11 + det*map.g01*map.g11;
  const real G3_1_0 = det*map.g10*map.g00 + det*map.g10*map.g00 + det*map.g11*map.g01 + det*map.g11*map.g01 + det*map.g11*map.g01 + det*map.g11*map.g01;
  const real G3_1_1 = det*map.g10*map.g10 + det*map.g10*map.g10 + det*map.g11*map.g11 + det*map.g11*map.g11 + det*map.g11*map.g11 + det*map.g11*map.g11;

  // Compute element tensor
  block[0] = 1.249999999999999e-01*G0_0_0 + 1.249999999999999e-01*G0_0_1 + 1.249999999999999e-01*G0_1_0 + 1.249999999999999e-01*G0_1_1;
  block[1] = -1.249999999999999e-01*G0_0_0 - 1.249999999999999e-01*G0_1_0;
  block[2] = -1.249999999999999e-01*G0_0_1 - 1.249999999999999e-01*G0_1_1;
  block[3] = 1.249999999999999e-01*G1_0_0 + 1.249999999999999e-01*G1_0_1 + 1.249999999999999e-01*G1_1_0 + 1.249999999999999e-01*G1_1_1;
  block[4] = -1.249999999999999e-01*G1_0_0 - 1.249999999999999e-01*G1_1_0;
  block[5] = -1.249999999999999e-01*G1_0_1 - 1.249999999999999e-01*G1_1_1;
  block[6] = -1.249999999999999e-01*G0_0_0 - 1.249999999999999e-01*G0_0_1;
  block[7] = 1.249999999999999e-01*G0_0_0;
  block[8] = 1.249999999999999e-01*G0_0_1;
  block[9] = -1.249999999999999e-01*G1_0_0 - 1.249999999999999e-01*G1_0_1;
  block[10] = 1.249999999999999e-01*G1_0_0;
  block[11] = 1.249999999999999e-01*G1_0_1;
  block[12] = -1.249999999999999e-01*G0_1_0 - 1.249999999999999e-01*G0_1_1;
  block[13] = 1.249999999999999e-01*G0_1_0;
  block[14] = 1.249999999999999e-01*G0_1_1;
  block[15] = -1.249999999999999e-01*G1_1_0 - 1.249999999999999e-01*G1_1_1;
  block[16] = 1.249999999999999e-01*G1_1_0;
  block[17] = 1.249999999999999e-01*G1_1_1;
  block[18] = 1.249999999999999e-01*G2_0_0 + 1.249999999999999e-01*G2_0_1 + 1.249999999999999e-01*G2_1_0 + 1.249999999999999e-01*G2_1_1;
  block[19] = -1.249999999999999e-01*G2_0_0 - 1.249999999999999e-01*G2_1_0;
  block[20] = -1.249999999999999e-01*G2_0_1 - 1.249999999999999e-01*G2_1_1;
  block[21] = 1.249999999999999e-01*G3_0_0 + 1.249999999999999e-01*G3_0_1 + 1.249999999999999e-01*G3_1_0 + 1.249999999999999e-01*G3_1_1;
  block[22] = -1.249999999999999e-01*G3_0_0 - 1.249999999999999e-01*G3_1_0;
  block[23] = -1.249999999999999e-01*G3_0_1 - 1.249999999999999e-01*G3_1_1;
  block[24] = -1.249999999999999e-01*G2_0_0 - 1.249999999999999e-01*G2_0_1;
  block[25] = 1.249999999999999e-01*G2_0_0;
  block[26] = 1.249999999999999e-01*G2_0_1;
  block[27] = -1.249999999999999e-01*G3_0_0 - 1.249999999999999e-01*G3_0_1;
  block[28] = 1.249999999999999e-01*G3_0_0;
  block[29] = 1.249999999999999e-01*G3_0_1;
  block[30] = -1.249999999999999e-01*G2_1_0 - 1.249999999999999e-01*G2_1_1;
  block[31] = 1.249999999999999e-01*G2_1_0;
  block[32] = 1.249999999999999e-01*G2_1_1;
  block[33] = -1.249999999999999e-01*G3_1_0 - 1.249999999999999e-01*G3_1_1;
  block[34] = 1.249999999999999e-01*G3_1_0;
  block[35] = 1.249999999999999e-01*G3_1_1;
}

// No contribution from the boundary
bool BilinearForm::boundary_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map, real det, unsigned int facet) const {}

// No contribution from interior boundaries
bool BilinearForm::interior_boundary_contribution() const { return false; }

void BilinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const {}

} }

#endif
