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
  parameters["tolerance"] = 5.0;

  // Access values
  int maxiter = parameters["max iterations"];
  double tol = parameters["tolerance"];
  cout << "maxiter = " << maxiter << endl;
  cout << "tol = " << tol << endl;

  // Print info
  parameters.info();

  return 0;
} 
