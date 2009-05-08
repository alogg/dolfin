#include <dolfin.h>
#include <dolfin/parameter/Parameters.h>

using namespace dolfin;

int main (int argc, char* argv[])
{
  // Application parameter database
  Parameters application_parameters("application parameters");

  // Set application parameters
  application_parameters.add("foo", 1.0);
  application_parameters.add("bar", 100);

  // Solver parameter database
  Parameters solver_parameters("solver parameters");

  // Set solver parameters
  solver_parameters.add("max iterations", 100);
  solver_parameters.add("tolerance", 1e-16);
  solver_parameters.add("relative tolerance", 1e-16, 1e-16, 1.0);
  
  // Set range
  solver_parameters("max iterations").set_range(0, 1000);  

  // Set values
  solver_parameters("max iterations") = 500;
  solver_parameters("relative tolerance") = 0.1;

  // Set solver parameters as nested parameters of application parameters
  application_parameters.add(solver_parameters);

  // Access values
  double foo = application_parameters("foo");
  int bar = application_parameters("bar");
  double tol = application_parameters["solver parameters"]("tolerance");
  
  // Silly hack to prevent warning from GCC about unused variables
  foo += 1; bar += 1; tol += 1;

  // Print parameters
  application_parameters.print();

  return 0;
} 
