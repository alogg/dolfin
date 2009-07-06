from dolfin import *

# Application parameter database
app_params = Parameters("application_parameters",
                           foo=1.0,
                           bar=100,
                           pc="amg",
                           solver = Parameters("solver",
                                                  max_iterations = 100,
                                                  tolerance = 1e-16,
                                                  relative_tolerance = (1e-16, 1e-16, 1.0),
                                                  prec = ("ilu",["ilu","amg","icc","sor"])
                                                  )
                           )
# Set range
app_params.solver.set_range("max_iterations",0, 1000)  

# Set values
app_params.solver.max_iterations = 500
app_params.solver.relative_tolerance = 0.1

# Parse command-line options
app_params.parse()

# Access values
foo = app_params.foo
bar = app_params["bar"]
tol = app_params.solver["tolerance"]

# Print parameters
logger.info(app_params)

# Test parameters for Krylov solver
solver = KrylovSolver()
solver.parameters["relative_tolerance"] = 1e-20
logger.info("")
logger.info(solver.parameters)


# Test a the use of update

# Update using an other Parameters
subset1 = Parameters("subset",
                       foo=3.0,
                       nested_subset = Parameters("solver",
                                                     max_iterations=850))

app_params.update(subset1)
logger.info("")
logger.info(app_params)

# Update using a dict
subset2 = dict(foo =1.5,solver = dict(max_iterations=50))
app_params.update(subset2)

logger.info("")
logger.info(app_params)

print "\nOption string representation of app_params:"
print app_params.option_string()
