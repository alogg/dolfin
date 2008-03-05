#!/bin/sh

# Comparing Metis and Scotch with varying number of cells and partitions
echo "Running benchmark 2"
echo "Running time is 1 minutes 10 seconds on a Sempron 3300+"
./benchmark metis1 unitcube4.xml 4 100000 metis_cube_4part.dat
./benchmark metis1 unitcube4.xml 16 100000 metis_cube_16part.dat
./benchmark metis1 unitcube4.xml 64 100000 metis_cube_64part.dat

./benchmark scotch unitcube4.xml 4 100000 scotch_cube_4part.dat
./benchmark scotch unitcube4.xml 16 100000 scotch_cube_16part.dat
./benchmark scotch unitcube4.xml 64 100000 scotch_cube_64part.dat
