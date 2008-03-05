#!/bin/sh

# Comparing Metis and Scotch with varying number of cells and partitions
echo "Running benchmark 1"
echo "Running time is 1 minute 14 seconds on a Sempron 3300+"
./benchmark metis1 unitsquare40.xml 4 600000 metis_square_4part.dat
./benchmark metis1 unitsquare40.xml 16 600000 metis_square_16part.dat
./benchmark metis1 unitsquare40.xml 64 600000 metis_square_64part.dat

./benchmark scotch unitsquare40.xml 4 600000 scotch_square_4part.dat
./benchmark scotch unitsquare40.xml 16 600000 scotch_square_16part.dat
./benchmark scotch unitsquare40.xml 64 600000 scotch_square_64part.dat
