#include <dolfin.h>

using namespace dolfin;

void test_parameters(int argc, char* argv[])
{
  // Application parameter database
  Parameters application_parameters("application_parameters");

  // Set application parameters
  application_parameters.add("foo", 1.0);
  application_parameters.add("bar", 100);
  application_parameters.add("pc", "amg");

  // Solver parameter database
  Parameters solver_parameters("solver_parameters");

  // Set solver parameters
  solver_parameters.add("max_iterations", 100);
  solver_parameters.add("tolerance", 1e-16);
  solver_parameters.add("relative_tolerance", 1e-16, 1e-16, 1.0);

  // Set range
  solver_parameters("max_iterations").set_range(0, 1000);

  // Set values
  solver_parameters("max_iterations") = 500;
  solver_parameters("relative_tolerance") = 0.1;

  // Set solver parameters as nested parameters of application parameters
  application_parameters.add(solver_parameters);

  // Parse command-line options
  application_parameters.parse(argc, argv);

  // Access values
  double foo = application_parameters("foo");
  int bar = application_parameters("bar");
  double tol = application_parameters["solver_parameters"]("tolerance");

  // Silly hack to prevent warning from GCC about unused variables
  foo += 1; bar += 1; tol += 1;

  // Print parameters
  info(application_parameters);

  // Test parameters for Krylov solver
  KrylovSolver solver;
  solver.parameters("relative_tolerance") = 1e-20;
  info("");
  info(solver.parameters);

  // Solver parameter database to be used together with update
  Parameters parameter_subset("parameter_subset");
  parameter_subset.add("foo", 3.0);

  Parameters nested_subset("solver_parameters");
  nested_subset.add("max_iterations", 850);

  parameter_subset.add(nested_subset);

  application_parameters.update(parameter_subset);

  // Print parameters
  info("");
  info(parameter_subset);
  info("");
  info(application_parameters);
}

int main (int argc, char* argv[])
{
  UnitSquare mesh(3, 3);
  mesh.init(1, 2);

  dolfin::cout << mesh << dolfin::endl;

  //info(mesh, true);

  return 0;
}
