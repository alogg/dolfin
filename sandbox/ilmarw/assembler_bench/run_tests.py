import os

from dolfin import *

assembly_tests = ['Elasticity_3D',
                  'ICNS_3D_Momentum',
                  'Laplace_2D',
                  'Stokes_2D_Stab',
                  'Stokes_2D_TH']

sizes = [32,
         32,
         256,
         256,
         256]

#print dolfin_set("linear algebra backend", "PETSc")

for i in range(0,len(assembly_tests),1):
    
    executable = './' + assembly_tests[i] + '/assembler_test'
    if os.path.isfile(executable):
        cmd = executable + ' ' + str(sizes[i])
        os.system(cmd)
    else:
        print 'File ' + executable + ' not found (did you scons?)'


