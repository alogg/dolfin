# Run this with mpirun -n 4 test_barrier.py

from dolfin import *

from time import sleep

# Sleeping for a second

if MPI.process_number() == 0:
    print "Process 0 sleeping for 5 seconds"
    sleep(5)
else:
    print "Process %d waiting for process 0 to finish" % MPI.process_number()

MPI.barrier()

print "Process %d done" % MPI.process_number()
