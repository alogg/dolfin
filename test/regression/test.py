"""Run all demos."""

__author__ = "Ilmar Wilbers (ilmarw@simula.no)"
__date__ = "2008-04-08 -- 2008-08-17"
__copyright__ = "Copyright (C) 2008 Ilmar Wilbers"
__license__  = "GNU LGPL Version 2.1"

# Modified by Anders Logg, 2008-2009.
# Modified by Johannes Ring, 2009.

import sys, os, re
import platform
from time import time
from subprocess import Popen, PIPE, STDOUT
from dolfin.utils import getstatusoutput

# Location of all demos
demodir = os.path.join(os.curdir, "..", "..", "demo")

# Demos to run
cppdemos = []
pydemos = []
for dpath, dnames, fnames in os.walk(demodir):
    if os.path.basename(dpath) == 'cpp':
        if os.path.isfile(os.path.join(dpath, 'SConstruct')):
            cppdemos.append(dpath)
    elif os.path.basename(dpath) == 'python':
        if os.path.isfile(os.path.join(dpath, 'demo.py')):
            pydemos.append(dpath)

# Set non-interactive
os.putenv('DOLFIN_NOPLOT', '1')

print "Running all demos (non-interactively)"
print ""
print "Found %d C++ demos" % len(cppdemos)
print "Found %d Python demos" % len(pydemos)
print ""
import pprint

# Remove demos that are known not to work (FIXME's)
pydemos.remove(os.path.join(demodir, 'ode', 'aliev-panfilov', 'python'))
pydemos.remove(os.path.join(demodir, 'ode', 'lorenz', 'python'))

# Push slow demos to the end
pyslow = []
cppslow = []
for s in pyslow:
    pydemos.remove(s)
    pydemos.append(s)
for s in cppslow:
    cppdemos.remove(s)
    cppdemos.append(s)

# Remove overly slow demos
pydemos.remove(os.path.join(demodir, 'pde', 'cahn-hilliard', 'python'))
cppdemos.remove(os.path.join(demodir, 'pde', 'cahn-hilliard', 'cpp'))

# Demos that need command-line arguments are treated separately
pydemos.remove(os.path.join(demodir, 'quadrature', 'python'))
cppdemos.remove(os.path.join(demodir, 'quadrature', 'cpp'))
cppdemos.remove(os.path.join(demodir, 'ode', 'method-weights', 'cpp'))
cppdemos.remove(os.path.join(demodir, 'ode', 'stiff', 'cpp'))

failed = []
timing = []

# Check if we should run only Python tests, use for quick testing
if len(sys.argv) == 2 and sys.argv[1] == "--only-python":
    only_python = True
else:
    only_python = False

# Check if we should skip C++ demos
if only_python:
    print "Skipping C++ demos"
    cppdemos = []

# Run in serial, then in parallel
for prefix in ["", "mpirun -n 2 "]:

    # Run C++ demos
    for demo in cppdemos:
        print "----------------------------------------------------------------------"
        print "Running C++ demo %s%s" % (prefix, demo)
        print ""
        cppdemo_ext = ''
        if platform.system() == 'Windows':
            cppdemo_ext = '.exe'
        if os.path.isfile(os.path.join(demo, 'demo' + cppdemo_ext)):
            t1 = time()
            output = getstatusoutput("cd %s && %s .%sdemo%s" % \
                                         (demo, prefix, os.path.sep, cppdemo_ext))
            t2 = time()
            timing += [(t2 - t1, demo)]
            if output[0] == 0:
                print "OK"
            elif output[0] == 10: # Failing but exiting gracefully
                print "ok (graceful exit on fail)"
            else:
                print "*** Failed"
                failed += [(demo, "C++", prefix, output[1])]
        else:
            print "*** Warning: missing demo"

    # Run Python demos
    for demo in pydemos:
        print "----------------------------------------------------------------------"
        print "Running Python demo %s%s" % (prefix, demo)
        print ""
        if os.path.isfile(os.path.join(demo, 'demo.py')):
            t1 = time()
            output = getstatusoutput("cd %s && %s python demo.py" % (demo, prefix))
            t2 = time()
            timing += [(t2 - t1, demo)]
            if output[0] == 0:
                print "OK"
            elif output[0] == 10: # Failing but exiting gracefully
                print "ok (graceful exit on fail)"
            else:
                print "*** Failed"
                failed += [(demo, "Python", prefix, output[1])]
        else:
            print "*** Warning: missing demo"

# Print summary of time to run demos
timing.sort()
print ""
print "Time to run demos:"
print "\n".join("%.2fs: %s" % t for t in timing)

total_no_demos = len(pydemos)
if not only_python:
    total_no_demos += len(cppdemos)

# Print output for failed tests
print ""
if len(failed) > 0:
    print "%d demo(s) out of %d failed, see demo.log for details." % \
          (len(failed), total_no_demos)
    file = open("demo.log", "w")
    for (test, interface, prefix, output) in failed:
        file.write("----------------------------------------------------------------------\n")
        file.write("%s%s (%s)\n" % (prefix, test, interface))
        file.write("\n")
        file.write(output)
        file.write("\n")
        file.write("\n")
else:
    print "All demos checked: OK"

# Return error code if tests failed
sys.exit(len(failed) != 0)
