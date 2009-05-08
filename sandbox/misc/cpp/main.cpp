#include <dolfin.h>
#include <dolfin/parameter/Parameters.h>

using namespace dolfin;

int main (int argc, char* argv[])
{
  // Create parameter database
  Parameters parameters("my parameters");

  // Define parameters
  parameters.add("max iterations", 100);
  parameters.add("tolerance", 1e-16);
  parameters.add("relative tolerance", 1e-16, 1e-16, 1.0);

  // Set range
  parameters["max iterations"].set_range(0, 1000);

  // Set values
  parameters["max iterations"] = 500;
  parameters["relative tolerance"] = 0.1;

  // Access values
  int maxiter = parameters["max iterations"];
  double tol = parameters["tolerance"];
  
  // Silly hack to prevent warning from GCC about unused variables
  maxiter++; tol += 1.0;

  // Print parameters
  parameters.print();

  return 0;
} 
