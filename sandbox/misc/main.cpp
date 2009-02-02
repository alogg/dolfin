// Place for random tests

#include <dolfin.h>

using namespace dolfin;

int main()
{
  // Create meshes
  //UnitCircle mesh0(10);
  UnitSquare mesh1(5,5);  
  UnitSquare mesh0(2,2);  
    
  cout << "***************HEJ**************" << endl;

  // Move and scale circle
   MeshGeometry& g = mesh0.geometry();
  for (VertexIterator vertex(mesh0); !vertex.end(); ++vertex)
  {
//     double* x = g.x(vertex->index());
//     x[0] = 0.5*x[0] + 1.05;
//     x[1] = 0.5*x[1] + 1.05;
  
    double* x = g.x(vertex->index());
    x[0] = x[0] + 0.75;
    x[1] = x[1] + 0.55; //0.60
  }

  // Compute intersection with boundary of square
  BoundaryMesh boundary1(mesh1);
  Array<unsigned int> cells0;
  mesh0.intersection(boundary1, cells0);

  // File f0("mesh0.pvd");
  //File f1("mesh1.pvd");
  //f0 << mesh0;
  //f1 << mesh1;
}


