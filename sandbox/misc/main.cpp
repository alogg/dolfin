#include <dolfin.h>

using namespace dolfin;

int main (int argc, char **argv) {

   UnitCube mesh(10, 10, 10);

   FiniteElement fe("FiniteElement('Lagrange', 'tetrahedron', 1)");
   DofMap dofmap("FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)", mesh);

   FunctionSpace V(mesh, fe, dofmap);
   Function j( V );

   j.vector() = 1.0;

   Point point(0.5, 0.5, 0.5);
   double x[3];

   j.eval(x, point.coordinates());
   std::cout << x[0] << std::endl;

   return 0;
} 
