// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.5.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __FSIDOTSIGMA3D_H
#define __FSIDOTSIGMA3D_H

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

namespace dolfin { namespace FSIDotSigma3D {

/// This class contains the form to be evaluated, including
/// contributions from the interior and boundary of the domain.

class BilinearForm : public dolfin::BilinearForm
{
public:

  class TestElement;

  class TrialElement;

  BilinearForm(const real& c0);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map, real det) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, real det, unsigned int facet) const;

  bool interior_boundary_contribution() const;

  void eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const;

private:

  const real& c0;

};

class BilinearForm::TestElement : public dolfin::FiniteElement
{
public:

  TestElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    tensordims = new unsigned int [1];
    tensordims[0] = 9;

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
    return 9;
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
    nodes[0] = cell.index();
    int offset = mesh.topology().size(3);
    nodes[1] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[2] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[3] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[4] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[5] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[6] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[7] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[8] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[1] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[2] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[3] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[4] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[5] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[6] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[7] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[8] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    components[0] = 0;
    components[1] = 1;
    components[2] = 2;
    components[3] = 3;
    components[4] = 4;
    components[5] = 5;
    components[6] = 6;
    components[7] = 7;
    components[8] = 8;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(3);
    vertex_nodes[1] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[2] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[3] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[4] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[5] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[6] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[7] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[8] = offset + vertex;
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
    FiniteElementSpec s("Discontinuous vector Lagrange", "tetrahedron", 0, 9);
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
    tensordims[0] = 9;

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
    return 9;
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
    nodes[0] = cell.index();
    int offset = mesh.topology().size(3);
    nodes[1] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[2] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[3] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[4] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[5] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[6] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[7] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[8] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[1] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[2] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[3] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[4] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[5] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[6] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[7] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[8] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    components[0] = 0;
    components[1] = 1;
    components[2] = 2;
    components[3] = 3;
    components[4] = 4;
    components[5] = 5;
    components[6] = 6;
    components[7] = 7;
    components[8] = 8;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(3);
    vertex_nodes[1] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[2] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[3] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[4] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[5] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[6] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[7] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[8] = offset + vertex;
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
    FiniteElementSpec s("Discontinuous vector Lagrange", "tetrahedron", 0, 9);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

BilinearForm::BilinearForm(const real& c0) : dolfin::BilinearForm(0), c0(c0)
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
  const real G0_ = det*c0;
  const real G1_ = det*c0;
  const real G2_ = det*c0;
  const real G3_ = det*c0;
  const real G4_ = det*c0;
  const real G5_ = det*c0;
  const real G6_ = det*c0;
  const real G7_ = det*c0;
  const real G8_ = det*c0;

  // Compute element tensor
  block[0] = 1.666666666666665e-01*G0_;
  block[1] = 0.000000000000000e+00;
  block[2] = 0.000000000000000e+00;
  block[3] = 0.000000000000000e+00;
  block[4] = 0.000000000000000e+00;
  block[5] = 0.000000000000000e+00;
  block[6] = 0.000000000000000e+00;
  block[7] = 0.000000000000000e+00;
  block[8] = 0.000000000000000e+00;
  block[9] = 0.000000000000000e+00;
  block[10] = 1.666666666666665e-01*G3_;
  block[11] = 0.000000000000000e+00;
  block[12] = 0.000000000000000e+00;
  block[13] = 0.000000000000000e+00;
  block[14] = 0.000000000000000e+00;
  block[15] = 0.000000000000000e+00;
  block[16] = 0.000000000000000e+00;
  block[17] = 0.000000000000000e+00;
  block[18] = 0.000000000000000e+00;
  block[19] = 0.000000000000000e+00;
  block[20] = 1.666666666666665e-01*G6_;
  block[21] = 0.000000000000000e+00;
  block[22] = 0.000000000000000e+00;
  block[23] = 0.000000000000000e+00;
  block[24] = 0.000000000000000e+00;
  block[25] = 0.000000000000000e+00;
  block[26] = 0.000000000000000e+00;
  block[27] = 0.000000000000000e+00;
  block[28] = 0.000000000000000e+00;
  block[29] = 0.000000000000000e+00;
  block[30] = 1.666666666666665e-01*G1_;
  block[31] = 0.000000000000000e+00;
  block[32] = 0.000000000000000e+00;
  block[33] = 0.000000000000000e+00;
  block[34] = 0.000000000000000e+00;
  block[35] = 0.000000000000000e+00;
  block[36] = 0.000000000000000e+00;
  block[37] = 0.000000000000000e+00;
  block[38] = 0.000000000000000e+00;
  block[39] = 0.000000000000000e+00;
  block[40] = 1.666666666666665e-01*G4_;
  block[41] = 0.000000000000000e+00;
  block[42] = 0.000000000000000e+00;
  block[43] = 0.000000000000000e+00;
  block[44] = 0.000000000000000e+00;
  block[45] = 0.000000000000000e+00;
  block[46] = 0.000000000000000e+00;
  block[47] = 0.000000000000000e+00;
  block[48] = 0.000000000000000e+00;
  block[49] = 0.000000000000000e+00;
  block[50] = 1.666666666666665e-01*G7_;
  block[51] = 0.000000000000000e+00;
  block[52] = 0.000000000000000e+00;
  block[53] = 0.000000000000000e+00;
  block[54] = 0.000000000000000e+00;
  block[55] = 0.000000000000000e+00;
  block[56] = 0.000000000000000e+00;
  block[57] = 0.000000000000000e+00;
  block[58] = 0.000000000000000e+00;
  block[59] = 0.000000000000000e+00;
  block[60] = 1.666666666666665e-01*G2_;
  block[61] = 0.000000000000000e+00;
  block[62] = 0.000000000000000e+00;
  block[63] = 0.000000000000000e+00;
  block[64] = 0.000000000000000e+00;
  block[65] = 0.000000000000000e+00;
  block[66] = 0.000000000000000e+00;
  block[67] = 0.000000000000000e+00;
  block[68] = 0.000000000000000e+00;
  block[69] = 0.000000000000000e+00;
  block[70] = 1.666666666666665e-01*G5_;
  block[71] = 0.000000000000000e+00;
  block[72] = 0.000000000000000e+00;
  block[73] = 0.000000000000000e+00;
  block[74] = 0.000000000000000e+00;
  block[75] = 0.000000000000000e+00;
  block[76] = 0.000000000000000e+00;
  block[77] = 0.000000000000000e+00;
  block[78] = 0.000000000000000e+00;
  block[79] = 0.000000000000000e+00;
  block[80] = 1.666666666666665e-01*G8_;
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

  LinearForm(Function& w0, const real& c0);
  

  bool interior_contribution() const;

  void eval(real block[], const AffineMap& map, real det) const;

  bool boundary_contribution() const;

  void eval(real block[], const AffineMap& map, real det, unsigned int facet) const;

  bool interior_boundary_contribution() const;

  void eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const;

private:

  const real& c0;

};

class LinearForm::TestElement : public dolfin::FiniteElement
{
public:

  TestElement() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    tensordims = new unsigned int [1];
    tensordims[0] = 9;

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
    return 9;
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
    nodes[0] = cell.index();
    int offset = mesh.topology().size(3);
    nodes[1] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[2] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[3] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[4] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[5] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[6] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[7] = offset + cell.index();
    offset = offset + mesh.topology().size(3);
    nodes[8] = offset + cell.index();
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[1] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[2] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[3] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[4] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[5] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[6] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[7] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    points[8] = map(2.500000000000000e-01, 2.500000000000000e-01, 2.500000000000000e-01);
    components[0] = 0;
    components[1] = 1;
    components[2] = 2;
    components[3] = 3;
    components[4] = 4;
    components[5] = 5;
    components[6] = 6;
    components[7] = 7;
    components[8] = 8;
  }

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(3);
    vertex_nodes[1] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[2] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[3] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[4] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[5] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[6] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[7] = offset + vertex;
    offset = offset + mesh.topology().size(3);
    vertex_nodes[8] = offset + vertex;
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
    FiniteElementSpec s("Discontinuous vector Lagrange", "tetrahedron", 0, 9);
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
    tensordims = new unsigned int [1];
    tensordims[0] = 3;

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
    return 3;
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
    nodes[3] = cell.entities(0)[3];
    int offset = mesh.topology().size(0);
    nodes[4] = offset + cell.entities(0)[0];
    nodes[5] = offset + cell.entities(0)[1];
    nodes[6] = offset + cell.entities(0)[2];
    nodes[7] = offset + cell.entities(0)[3];
    offset = offset + mesh.topology().size(0);
    nodes[8] = offset + cell.entities(0)[0];
    nodes[9] = offset + cell.entities(0)[1];
    nodes[10] = offset + cell.entities(0)[2];
    nodes[11] = offset + cell.entities(0)[3];
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

  void vertexeval(uint vertex_nodes[], unsigned int vertex, const Mesh& mesh) const
  {
    // FIXME: Temporary fix for Lagrange elements
    vertex_nodes[0] = vertex;
    int offset = mesh.topology().size(0);
    vertex_nodes[1] = offset + vertex;
    offset = offset + mesh.topology().size(0);
    vertex_nodes[2] = offset + vertex;
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
    FiniteElementSpec s("Vector Lagrange", "tetrahedron", 1, 3);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

LinearForm::LinearForm(Function& w0, const real& c0) : dolfin::LinearForm(1), c0(c0)
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
  const real c0_10 = c[0][10];
  const real c0_11 = c[0][11];

  // Compute geometry tensors
  const real G0_0_0 = det*c0*c0_0*map.g00 + det*c0*c0_0*map.g00;
  const real G0_0_1 = det*c0*c0_0*map.g10 + det*c0*c0_0*map.g10;
  const real G0_0_2 = det*c0*c0_0*map.g20 + det*c0*c0_0*map.g20;
  const real G0_1_0 = det*c0*c0_1*map.g00 + det*c0*c0_1*map.g00;
  const real G0_2_1 = det*c0*c0_2*map.g10 + det*c0*c0_2*map.g10;
  const real G0_3_2 = det*c0*c0_3*map.g20 + det*c0*c0_3*map.g20;
  const real G1_0_0 = det*c0*c0_0*map.g01;
  const real G1_0_1 = det*c0*c0_0*map.g11;
  const real G1_0_2 = det*c0*c0_0*map.g21;
  const real G1_1_0 = det*c0*c0_1*map.g01;
  const real G1_2_1 = det*c0*c0_2*map.g11;
  const real G1_3_2 = det*c0*c0_3*map.g21;
  const real G2_4_0 = det*c0*c0_4*map.g00;
  const real G2_4_1 = det*c0*c0_4*map.g10;
  const real G2_4_2 = det*c0*c0_4*map.g20;
  const real G2_5_0 = det*c0*c0_5*map.g00;
  const real G2_6_1 = det*c0*c0_6*map.g10;
  const real G2_7_2 = det*c0*c0_7*map.g20;
  const real G3_0_0 = det*c0*c0_0*map.g02;
  const real G3_0_1 = det*c0*c0_0*map.g12;
  const real G3_0_2 = det*c0*c0_0*map.g22;
  const real G3_1_0 = det*c0*c0_1*map.g02;
  const real G3_2_1 = det*c0*c0_2*map.g12;
  const real G3_3_2 = det*c0*c0_3*map.g22;
  const real G4_8_0 = det*c0*c0_8*map.g00;
  const real G4_8_1 = det*c0*c0_8*map.g10;
  const real G4_8_2 = det*c0*c0_8*map.g20;
  const real G4_9_0 = det*c0*c0_9*map.g00;
  const real G4_10_1 = det*c0*c0_10*map.g10;
  const real G4_11_2 = det*c0*c0_11*map.g20;
  const real G5_4_0 = det*c0*c0_4*map.g00;
  const real G5_4_1 = det*c0*c0_4*map.g10;
  const real G5_4_2 = det*c0*c0_4*map.g20;
  const real G5_5_0 = det*c0*c0_5*map.g00;
  const real G5_6_1 = det*c0*c0_6*map.g10;
  const real G5_7_2 = det*c0*c0_7*map.g20;
  const real G6_0_0 = det*c0*c0_0*map.g01;
  const real G6_0_1 = det*c0*c0_0*map.g11;
  const real G6_0_2 = det*c0*c0_0*map.g21;
  const real G6_1_0 = det*c0*c0_1*map.g01;
  const real G6_2_1 = det*c0*c0_2*map.g11;
  const real G6_3_2 = det*c0*c0_3*map.g21;
  const real G7_4_0 = det*c0*c0_4*map.g01 + det*c0*c0_4*map.g01;
  const real G7_4_1 = det*c0*c0_4*map.g11 + det*c0*c0_4*map.g11;
  const real G7_4_2 = det*c0*c0_4*map.g21 + det*c0*c0_4*map.g21;
  const real G7_5_0 = det*c0*c0_5*map.g01 + det*c0*c0_5*map.g01;
  const real G7_6_1 = det*c0*c0_6*map.g11 + det*c0*c0_6*map.g11;
  const real G7_7_2 = det*c0*c0_7*map.g21 + det*c0*c0_7*map.g21;
  const real G8_4_0 = det*c0*c0_4*map.g02;
  const real G8_4_1 = det*c0*c0_4*map.g12;
  const real G8_4_2 = det*c0*c0_4*map.g22;
  const real G8_5_0 = det*c0*c0_5*map.g02;
  const real G8_6_1 = det*c0*c0_6*map.g12;
  const real G8_7_2 = det*c0*c0_7*map.g22;
  const real G9_8_0 = det*c0*c0_8*map.g01;
  const real G9_8_1 = det*c0*c0_8*map.g11;
  const real G9_8_2 = det*c0*c0_8*map.g21;
  const real G9_9_0 = det*c0*c0_9*map.g01;
  const real G9_10_1 = det*c0*c0_10*map.g11;
  const real G9_11_2 = det*c0*c0_11*map.g21;
  const real G10_8_0 = det*c0*c0_8*map.g00;
  const real G10_8_1 = det*c0*c0_8*map.g10;
  const real G10_8_2 = det*c0*c0_8*map.g20;
  const real G10_9_0 = det*c0*c0_9*map.g00;
  const real G10_10_1 = det*c0*c0_10*map.g10;
  const real G10_11_2 = det*c0*c0_11*map.g20;
  const real G11_0_0 = det*c0*c0_0*map.g02;
  const real G11_0_1 = det*c0*c0_0*map.g12;
  const real G11_0_2 = det*c0*c0_0*map.g22;
  const real G11_1_0 = det*c0*c0_1*map.g02;
  const real G11_2_1 = det*c0*c0_2*map.g12;
  const real G11_3_2 = det*c0*c0_3*map.g22;
  const real G12_8_0 = det*c0*c0_8*map.g01;
  const real G12_8_1 = det*c0*c0_8*map.g11;
  const real G12_8_2 = det*c0*c0_8*map.g21;
  const real G12_9_0 = det*c0*c0_9*map.g01;
  const real G12_10_1 = det*c0*c0_10*map.g11;
  const real G12_11_2 = det*c0*c0_11*map.g21;
  const real G13_4_0 = det*c0*c0_4*map.g02;
  const real G13_4_1 = det*c0*c0_4*map.g12;
  const real G13_4_2 = det*c0*c0_4*map.g22;
  const real G13_5_0 = det*c0*c0_5*map.g02;
  const real G13_6_1 = det*c0*c0_6*map.g12;
  const real G13_7_2 = det*c0*c0_7*map.g22;
  const real G14_8_0 = det*c0*c0_8*map.g02 + det*c0*c0_8*map.g02;
  const real G14_8_1 = det*c0*c0_8*map.g12 + det*c0*c0_8*map.g12;
  const real G14_8_2 = det*c0*c0_8*map.g22 + det*c0*c0_8*map.g22;
  const real G14_9_0 = det*c0*c0_9*map.g02 + det*c0*c0_9*map.g02;
  const real G14_10_1 = det*c0*c0_10*map.g12 + det*c0*c0_10*map.g12;
  const real G14_11_2 = det*c0*c0_11*map.g22 + det*c0*c0_11*map.g22;

  // Compute element tensor
  block[0] = -1.666666666666665e-01*G0_0_0 - 1.666666666666665e-01*G0_0_1 - 1.666666666666664e-01*G0_0_2 + 1.666666666666665e-01*G0_1_0 + 1.666666666666665e-01*G0_2_1 + 1.666666666666665e-01*G0_3_2;
  block[1] = -1.666666666666665e-01*G5_4_0 - 1.666666666666665e-01*G5_4_1 - 1.666666666666664e-01*G5_4_2 + 1.666666666666665e-01*G5_5_0 + 1.666666666666665e-01*G5_6_1 + 1.666666666666665e-01*G5_7_2 - 1.666666666666665e-01*G6_0_0 - 1.666666666666665e-01*G6_0_1 - 1.666666666666664e-01*G6_0_2 + 1.666666666666665e-01*G6_1_0 + 1.666666666666665e-01*G6_2_1 + 1.666666666666665e-01*G6_3_2;
  block[2] = -1.666666666666665e-01*G10_8_0 - 1.666666666666665e-01*G10_8_1 - 1.666666666666664e-01*G10_8_2 + 1.666666666666665e-01*G10_9_0 + 1.666666666666665e-01*G10_10_1 + 1.666666666666665e-01*G10_11_2 - 1.666666666666665e-01*G11_0_0 - 1.666666666666665e-01*G11_0_1 - 1.666666666666664e-01*G11_0_2 + 1.666666666666665e-01*G11_1_0 + 1.666666666666665e-01*G11_2_1 + 1.666666666666665e-01*G11_3_2;
  block[3] = -1.666666666666665e-01*G1_0_0 - 1.666666666666665e-01*G1_0_1 - 1.666666666666664e-01*G1_0_2 + 1.666666666666665e-01*G1_1_0 + 1.666666666666665e-01*G1_2_1 + 1.666666666666665e-01*G1_3_2 - 1.666666666666665e-01*G2_4_0 - 1.666666666666665e-01*G2_4_1 - 1.666666666666664e-01*G2_4_2 + 1.666666666666665e-01*G2_5_0 + 1.666666666666665e-01*G2_6_1 + 1.666666666666665e-01*G2_7_2;
  block[4] = -1.666666666666665e-01*G7_4_0 - 1.666666666666665e-01*G7_4_1 - 1.666666666666664e-01*G7_4_2 + 1.666666666666665e-01*G7_5_0 + 1.666666666666665e-01*G7_6_1 + 1.666666666666665e-01*G7_7_2;
  block[5] = -1.666666666666665e-01*G12_8_0 - 1.666666666666665e-01*G12_8_1 - 1.666666666666664e-01*G12_8_2 + 1.666666666666665e-01*G12_9_0 + 1.666666666666665e-01*G12_10_1 + 1.666666666666665e-01*G12_11_2 - 1.666666666666665e-01*G13_4_0 - 1.666666666666665e-01*G13_4_1 - 1.666666666666664e-01*G13_4_2 + 1.666666666666665e-01*G13_5_0 + 1.666666666666665e-01*G13_6_1 + 1.666666666666665e-01*G13_7_2;
  block[6] = -1.666666666666665e-01*G3_0_0 - 1.666666666666665e-01*G3_0_1 - 1.666666666666664e-01*G3_0_2 + 1.666666666666665e-01*G3_1_0 + 1.666666666666665e-01*G3_2_1 + 1.666666666666665e-01*G3_3_2 - 1.666666666666665e-01*G4_8_0 - 1.666666666666665e-01*G4_8_1 - 1.666666666666664e-01*G4_8_2 + 1.666666666666665e-01*G4_9_0 + 1.666666666666665e-01*G4_10_1 + 1.666666666666665e-01*G4_11_2;
  block[7] = -1.666666666666665e-01*G8_4_0 - 1.666666666666665e-01*G8_4_1 - 1.666666666666664e-01*G8_4_2 + 1.666666666666665e-01*G8_5_0 + 1.666666666666665e-01*G8_6_1 + 1.666666666666665e-01*G8_7_2 - 1.666666666666665e-01*G9_8_0 - 1.666666666666665e-01*G9_8_1 - 1.666666666666664e-01*G9_8_2 + 1.666666666666665e-01*G9_9_0 + 1.666666666666665e-01*G9_10_1 + 1.666666666666665e-01*G9_11_2;
  block[8] = -1.666666666666665e-01*G14_8_0 - 1.666666666666665e-01*G14_8_1 - 1.666666666666664e-01*G14_8_2 + 1.666666666666665e-01*G14_9_0 + 1.666666666666665e-01*G14_10_1 + 1.666666666666665e-01*G14_11_2;
}

// No contribution from the boundary
bool LinearForm::boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map, real det, unsigned int facet) const {}

// No contribution from interior boundaries
bool LinearForm::interior_boundary_contribution() const { return false; }

void LinearForm::eval(real block[], const AffineMap& map0, const AffineMap& map1, real det, unsigned int facet0, unsigned int facet1, unsigned int alignment) const {}

} }

#endif
