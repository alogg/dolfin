echo "Poisson on UnitCube(30, 30, 30)" > poisson_unitcube30.dat
echo "Processes Time" >> poisson_unitcube30.dat
mpirun -n 1 ./dolfin-pmesh-test --cells 30 --num_iterations 1 --resultfile poisson_unitcube30.dat
mpirun -n 2 ./dolfin-pmesh-test --cells 30 --num_iterations 1 --resultfile poisson_unitcube30.dat
mpirun -n 3 ./dolfin-pmesh-test --cells 30 --num_iterations 1 --resultfile poisson_unitcube30.dat
mpirun -n 4 ./dolfin-pmesh-test --cells 30 --num_iterations 1 --resultfile poisson_unitcube30.dat
