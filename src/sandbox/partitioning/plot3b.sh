#!/bin/sh

if [ -r metis_torso_4part.dat ] 
then
python plot.py Torso edgecut metis_torso_4part.dat metis_torso_16part.dat metis_torso_64part.dat scotch_torso_4part.dat scotch_torso_16part.dat scotch_torso_64part.dat
exit
fi

echo "Run bench3.sh before running this"
