// Place for random tests

#include <dolfin.h>

using namespace dolfin;

int main()
{
  Table table;

  table("uBLAS",  "Assemble") = 0.010;
  table("uBLAS",  "Solve")    = 0.020;
  table("PETSc",  "Assemble") = 0.011;
  table("PETSc",  "Solve")    = 0.019;
  table("Epetra", "Assemble") = 0.012;
  table("Epetra", "Solve")    = 0.018;

  table.disp();
}
