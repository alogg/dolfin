#!/bin/sh

# Comparing Metis and Scotch with varying number of cells and partitions
echo "Running benchmark 1"
echo "Running time is 1 minute 14 seconds on a Sempron 3300+"
./demo metis1 unitsquare40.xml 4 600000 metis_square_4part.dat
./demo metis1 unitsquare40.xml 16 600000 metis_square_16part.dat
./demo metis1 unitsquare40.xml 64 600000 metis_square_64part.dat

./demo scotch unitsquare40.xml 4 600000 scotch_square_4part.dat
./demo scotch unitsquare40.xml 16 600000 scotch_square_16part.dat
./demo scotch unitsquare40.xml 64 600000 scotch_square_64part.dat
