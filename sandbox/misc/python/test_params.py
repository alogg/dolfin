from dolfin import *
from dolfin.cpp import NewParameters
import sys

# Application parameter database
application_parameters = NewParameters("application_parameters")

# Set application parameters
application_parameters.add("foo", 1.0);
application_parameters.add("bar", 100)
application_parameters.add("pc", "amg")

# Solver parameter database
solver_parameters = NewParameters ("solver_parameters")

# Set solver parameters
solver_parameters.add("max_iterations", 100)
solver_parameters.add("tolerance", 1e-16)
solver_parameters.add("relative_tolerance", 1e-16, 1e-16, 1.0)
solver_parameters.add("prec","ilu",["ilu","amg","icc","sor"])
  
# Set range
solver_parameters.set_range("max_iterations",0, 1000)  

# Set values
solver_parameters["max_iterations"] = 500
solver_parameters["relative_tolerance"] = 0.1

# Set solver parameters as nested parameters of application parameters
application_parameters.add(solver_parameters)

# Parse command-line options
application_parameters._parse(sys.argv)

# Access values
foo = application_parameters["foo"]
bar = application_parameters["bar"]
tol = application_parameters["solver_parameters"]["tolerance"]

# Print parameters
logger.info(application_parameters)

# Test parameters for Krylov solver
solver = KrylovSolver()
solver.parameters["relative_tolerance"] = 1e-20
logger.info("")
logger.info(solver.parameters)
