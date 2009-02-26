#include <dolfin.h>

using namespace dolfin;

int main (int argc, char **argv) {

   UnitCube mesh(10,10,10);

   FiniteElement fe("FiniteElement('Lagrange', 'tetrahedron', 1)");
   DofMap dofmap("FFC dof map for FiniteElement('Lagrange', 'tetrahedron', 1)", mesh);

   Function j( FunctionSpace(mesh, fe, dofmap) );

   j.vector() = 1.0;


   Point point( 0.5, 0.5, 0.5);
   Point value;

   j.eval((double *)value.coordinates(), point.coordinates());

   std::cout << value[0] << std::endl;

   return 0;
} 
