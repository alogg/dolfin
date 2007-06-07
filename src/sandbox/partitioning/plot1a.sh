#!/bin/sh

if [ -r metis_square_4part.dat ] 
then
python plot.py UnitSquare time metis_square_4part.dat metis_square_16part.dat metis_square_64part.dat scotch_square_4part.dat scotch_square_16part.dat scotch_square_64part.dat
exit
fi

echo "Run bench1.sh before running this"
