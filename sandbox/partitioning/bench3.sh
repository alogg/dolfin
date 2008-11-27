#!/bin/sh

# Comparing Metis and Scotch with varying number of cells and partitions
echo "Running benchmark 3"
echo "Running time is 2 minutes 13 seconds on a Sempron 3300+"
./demo metis1 ../../data/meshes/torso.xml.gz 4 100000 metis_torso_4part.dat
./demo metis1 ../../data/meshes/torso.xml.gz 16 100000 metis_torso_16part.dat
./demo metis1 ../../data/meshes/torso.xml.gz 64 100000 metis_torso_64part.dat

./demo scotch ../../data/meshes/torso.xml.gz 4 100000 scotch_torso_4part.dat
./demo scotch ../../data/meshes/torso.xml.gz 16 100000 scotch_torso_16part.dat
./demo scotch ../../data/meshes/torso.xml.gz 64 100000 scotch_torso_64part.dat
