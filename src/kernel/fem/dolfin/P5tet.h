// Automatically generated by FFC, the FEniCS Form Compiler, version 0.3.3-dev.
// For further information, go to http://www/fenics.org/ffc/.
// Licensed under the GNU GPL Version 2.

#ifndef __P5TET_H
#define __P5TET_H

#include <dolfin/Mesh.h>
#include <dolfin/Cell.h>
#include <dolfin/Point.h>
#include <dolfin/AffineMap.h>
#include <dolfin/FiniteElement.h>
#include <dolfin/FiniteElementSpec.h>

namespace dolfin
{

class P5tet : public dolfin::FiniteElement
{
public:

  P5tet() : dolfin::FiniteElement(), tensordims(0), subelements(0)
  {
    // Element is scalar, don't need to initialize tensordims

    // Element is simple, don't need to initialize subelements
  }

  ~P5tet()
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
    return 56;
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
    static unsigned int edge_reordering_0[2][4] = {{0, 1, 2, 3}, {3, 2, 1, 0}};
    static unsigned int face_reordering_0[6][6] = {{0, 1, 2, 3, 4, 5}, {0, 3, 5, 1, 4, 2}, {5, 3, 0, 4, 1, 2}, {2, 1, 0, 4, 3, 5}, {2, 4, 5, 1, 3, 0}, {5, 4, 2, 3, 1, 0}};
    nodes[0] = cell.entities(0)[0];
    nodes[1] = cell.entities(0)[1];
    nodes[2] = cell.entities(0)[2];
    nodes[3] = cell.entities(0)[3];
    int alignment = cell.alignment(1, 0);
    int offset = mesh.topology().size(0);
    nodes[4] = offset + 4*cell.entities(1)[0] + edge_reordering_0[alignment][0];
    nodes[5] = offset + 4*cell.entities(1)[0] + edge_reordering_0[alignment][1];
    nodes[6] = offset + 4*cell.entities(1)[0] + edge_reordering_0[alignment][2];
    nodes[7] = offset + 4*cell.entities(1)[0] + edge_reordering_0[alignment][3];
    alignment = cell.alignment(1, 1);
    nodes[8] = offset + 4*cell.entities(1)[1] + edge_reordering_0[alignment][0];
    nodes[9] = offset + 4*cell.entities(1)[1] + edge_reordering_0[alignment][1];
    nodes[10] = offset + 4*cell.entities(1)[1] + edge_reordering_0[alignment][2];
    nodes[11] = offset + 4*cell.entities(1)[1] + edge_reordering_0[alignment][3];
    alignment = cell.alignment(1, 2);
    nodes[12] = offset + 4*cell.entities(1)[2] + edge_reordering_0[alignment][0];
    nodes[13] = offset + 4*cell.entities(1)[2] + edge_reordering_0[alignment][1];
    nodes[14] = offset + 4*cell.entities(1)[2] + edge_reordering_0[alignment][2];
    nodes[15] = offset + 4*cell.entities(1)[2] + edge_reordering_0[alignment][3];
    alignment = cell.alignment(1, 3);
    nodes[16] = offset + 4*cell.entities(1)[3] + edge_reordering_0[alignment][0];
    nodes[17] = offset + 4*cell.entities(1)[3] + edge_reordering_0[alignment][1];
    nodes[18] = offset + 4*cell.entities(1)[3] + edge_reordering_0[alignment][2];
    nodes[19] = offset + 4*cell.entities(1)[3] + edge_reordering_0[alignment][3];
    alignment = cell.alignment(1, 4);
    nodes[20] = offset + 4*cell.entities(1)[4] + edge_reordering_0[alignment][0];
    nodes[21] = offset + 4*cell.entities(1)[4] + edge_reordering_0[alignment][1];
    nodes[22] = offset + 4*cell.entities(1)[4] + edge_reordering_0[alignment][2];
    nodes[23] = offset + 4*cell.entities(1)[4] + edge_reordering_0[alignment][3];
    alignment = cell.alignment(1, 5);
    nodes[24] = offset + 4*cell.entities(1)[5] + edge_reordering_0[alignment][0];
    nodes[25] = offset + 4*cell.entities(1)[5] + edge_reordering_0[alignment][1];
    nodes[26] = offset + 4*cell.entities(1)[5] + edge_reordering_0[alignment][2];
    nodes[27] = offset + 4*cell.entities(1)[5] + edge_reordering_0[alignment][3];
    alignment = cell.alignment(2, 0);
    offset = offset + 4*mesh.topology().size(1);
    nodes[28] = offset + 6*cell.entities(2)[0] + face_reordering_0[alignment][0];
    nodes[29] = offset + 6*cell.entities(2)[0] + face_reordering_0[alignment][1];
    nodes[30] = offset + 6*cell.entities(2)[0] + face_reordering_0[alignment][2];
    nodes[31] = offset + 6*cell.entities(2)[0] + face_reordering_0[alignment][3];
    nodes[32] = offset + 6*cell.entities(2)[0] + face_reordering_0[alignment][4];
    nodes[33] = offset + 6*cell.entities(2)[0] + face_reordering_0[alignment][5];
    alignment = cell.alignment(2, 1);
    nodes[34] = offset + 6*cell.entities(2)[1] + face_reordering_0[alignment][0];
    nodes[35] = offset + 6*cell.entities(2)[1] + face_reordering_0[alignment][1];
    nodes[36] = offset + 6*cell.entities(2)[1] + face_reordering_0[alignment][2];
    nodes[37] = offset + 6*cell.entities(2)[1] + face_reordering_0[alignment][3];
    nodes[38] = offset + 6*cell.entities(2)[1] + face_reordering_0[alignment][4];
    nodes[39] = offset + 6*cell.entities(2)[1] + face_reordering_0[alignment][5];
    alignment = cell.alignment(2, 2);
    nodes[40] = offset + 6*cell.entities(2)[2] + face_reordering_0[alignment][0];
    nodes[41] = offset + 6*cell.entities(2)[2] + face_reordering_0[alignment][1];
    nodes[42] = offset + 6*cell.entities(2)[2] + face_reordering_0[alignment][2];
    nodes[43] = offset + 6*cell.entities(2)[2] + face_reordering_0[alignment][3];
    nodes[44] = offset + 6*cell.entities(2)[2] + face_reordering_0[alignment][4];
    nodes[45] = offset + 6*cell.entities(2)[2] + face_reordering_0[alignment][5];
    alignment = cell.alignment(2, 3);
    nodes[46] = offset + 6*cell.entities(2)[3] + face_reordering_0[alignment][0];
    nodes[47] = offset + 6*cell.entities(2)[3] + face_reordering_0[alignment][1];
    nodes[48] = offset + 6*cell.entities(2)[3] + face_reordering_0[alignment][2];
    nodes[49] = offset + 6*cell.entities(2)[3] + face_reordering_0[alignment][3];
    nodes[50] = offset + 6*cell.entities(2)[3] + face_reordering_0[alignment][4];
    nodes[51] = offset + 6*cell.entities(2)[3] + face_reordering_0[alignment][5];
    offset = offset + 6*mesh.topology().size(2);
    nodes[52] = offset + 4*cell.index() + 0;
    nodes[53] = offset + 4*cell.index() + 1;
    nodes[54] = offset + 4*cell.index() + 2;
    nodes[55] = offset + 4*cell.index() + 3;
  }

  void pointmap(Point points[], unsigned int components[], const AffineMap& map) const
  {
    points[0] = map(0.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[1] = map(1.000000000000000e+00, 0.000000000000000e+00, 0.000000000000000e+00);
    points[2] = map(0.000000000000000e+00, 1.000000000000000e+00, 0.000000000000000e+00);
    points[3] = map(0.000000000000000e+00, 0.000000000000000e+00, 1.000000000000000e+00);
    points[4] = map(8.000000000000000e-01, 2.000000000000000e-01, 0.000000000000000e+00);
    points[5] = map(6.000000000000000e-01, 4.000000000000000e-01, 0.000000000000000e+00);
    points[6] = map(4.000000000000000e-01, 6.000000000000000e-01, 0.000000000000000e+00);
    points[7] = map(2.000000000000000e-01, 8.000000000000000e-01, 0.000000000000000e+00);
    points[8] = map(0.000000000000000e+00, 8.000000000000000e-01, 0.000000000000000e+00);
    points[9] = map(0.000000000000000e+00, 6.000000000000000e-01, 0.000000000000000e+00);
    points[10] = map(0.000000000000000e+00, 4.000000000000000e-01, 0.000000000000000e+00);
    points[11] = map(0.000000000000000e+00, 2.000000000000000e-01, 0.000000000000000e+00);
    points[12] = map(2.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
    points[13] = map(4.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
    points[14] = map(6.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
    points[15] = map(8.000000000000000e-01, 0.000000000000000e+00, 0.000000000000000e+00);
    points[16] = map(0.000000000000000e+00, 0.000000000000000e+00, 2.000000000000000e-01);
    points[17] = map(0.000000000000000e+00, 0.000000000000000e+00, 4.000000000000000e-01);
    points[18] = map(0.000000000000000e+00, 0.000000000000000e+00, 6.000000000000000e-01);
    points[19] = map(0.000000000000000e+00, 0.000000000000000e+00, 8.000000000000000e-01);
    points[20] = map(8.000000000000000e-01, 0.000000000000000e+00, 2.000000000000000e-01);
    points[21] = map(6.000000000000000e-01, 0.000000000000000e+00, 4.000000000000000e-01);
    points[22] = map(4.000000000000000e-01, 0.000000000000000e+00, 6.000000000000000e-01);
    points[23] = map(2.000000000000000e-01, 0.000000000000000e+00, 8.000000000000000e-01);
    points[24] = map(0.000000000000000e+00, 8.000000000000000e-01, 2.000000000000000e-01);
    points[25] = map(0.000000000000000e+00, 6.000000000000000e-01, 4.000000000000000e-01);
    points[26] = map(0.000000000000000e+00, 4.000000000000000e-01, 6.000000000000000e-01);
    points[27] = map(0.000000000000000e+00, 2.000000000000000e-01, 8.000000000000000e-01);
    points[28] = map(6.000000000000000e-01, 2.000000000000000e-01, 2.000000000000000e-01);
    points[29] = map(4.000000000000000e-01, 2.000000000000000e-01, 4.000000000000000e-01);
    points[30] = map(2.000000000000000e-01, 2.000000000000000e-01, 6.000000000000000e-01);
    points[31] = map(4.000000000000000e-01, 4.000000000000000e-01, 2.000000000000000e-01);
    points[32] = map(2.000000000000000e-01, 4.000000000000000e-01, 4.000000000000000e-01);
    points[33] = map(2.000000000000000e-01, 6.000000000000000e-01, 2.000000000000000e-01);
    points[34] = map(0.000000000000000e+00, 6.000000000000000e-01, 2.000000000000000e-01);
    points[35] = map(0.000000000000000e+00, 4.000000000000000e-01, 4.000000000000000e-01);
    points[36] = map(0.000000000000000e+00, 2.000000000000000e-01, 6.000000000000000e-01);
    points[37] = map(0.000000000000000e+00, 4.000000000000000e-01, 2.000000000000000e-01);
    points[38] = map(0.000000000000000e+00, 2.000000000000000e-01, 4.000000000000000e-01);
    points[39] = map(0.000000000000000e+00, 2.000000000000000e-01, 2.000000000000000e-01);
    points[40] = map(2.000000000000000e-01, 0.000000000000000e+00, 6.000000000000000e-01);
    points[41] = map(4.000000000000000e-01, 0.000000000000000e+00, 4.000000000000000e-01);
    points[42] = map(6.000000000000000e-01, 0.000000000000000e+00, 2.000000000000000e-01);
    points[43] = map(2.000000000000000e-01, 0.000000000000000e+00, 4.000000000000000e-01);
    points[44] = map(4.000000000000000e-01, 0.000000000000000e+00, 2.000000000000000e-01);
    points[45] = map(2.000000000000000e-01, 0.000000000000000e+00, 2.000000000000000e-01);
    points[46] = map(2.000000000000000e-01, 2.000000000000000e-01, 0.000000000000000e+00);
    points[47] = map(4.000000000000000e-01, 2.000000000000000e-01, 0.000000000000000e+00);
    points[48] = map(6.000000000000000e-01, 2.000000000000000e-01, 0.000000000000000e+00);
    points[49] = map(2.000000000000000e-01, 4.000000000000000e-01, 0.000000000000000e+00);
    points[50] = map(4.000000000000000e-01, 4.000000000000000e-01, 0.000000000000000e+00);
    points[51] = map(2.000000000000000e-01, 6.000000000000000e-01, 0.000000000000000e+00);
    points[52] = map(2.000000000000000e-01, 2.000000000000000e-01, 2.000000000000000e-01);
    points[53] = map(4.000000000000000e-01, 2.000000000000000e-01, 2.000000000000000e-01);
    points[54] = map(2.000000000000000e-01, 4.000000000000000e-01, 2.000000000000000e-01);
    points[55] = map(2.000000000000000e-01, 2.000000000000000e-01, 4.000000000000000e-01);
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
    components[10] = 0;
    components[11] = 0;
    components[12] = 0;
    components[13] = 0;
    components[14] = 0;
    components[15] = 0;
    components[16] = 0;
    components[17] = 0;
    components[18] = 0;
    components[19] = 0;
    components[20] = 0;
    components[21] = 0;
    components[22] = 0;
    components[23] = 0;
    components[24] = 0;
    components[25] = 0;
    components[26] = 0;
    components[27] = 0;
    components[28] = 0;
    components[29] = 0;
    components[30] = 0;
    components[31] = 0;
    components[32] = 0;
    components[33] = 0;
    components[34] = 0;
    components[35] = 0;
    components[36] = 0;
    components[37] = 0;
    components[38] = 0;
    components[39] = 0;
    components[40] = 0;
    components[41] = 0;
    components[42] = 0;
    components[43] = 0;
    components[44] = 0;
    components[45] = 0;
    components[46] = 0;
    components[47] = 0;
    components[48] = 0;
    components[49] = 0;
    components[50] = 0;
    components[51] = 0;
    components[52] = 0;
    components[53] = 0;
    components[54] = 0;
    components[55] = 0;
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
    FiniteElementSpec s("Lagrange", "tetrahedron", 5);
    return s;
  }
  
private:

  unsigned int* tensordims;
  FiniteElement** subelements;

};

}

#endif
