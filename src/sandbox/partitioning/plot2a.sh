#!/bin/sh

if [ -r metis_cube_4part.dat ] 
then
python plot.py UnitCube time metis_cube_4part.dat metis_cube_16part.dat metis_cube_64part.dat scotch_cube_4part.dat scotch_cube_16part.dat scotch_cube_64part.dat
exit
fi

echo "Run bench2.sh before running this"
