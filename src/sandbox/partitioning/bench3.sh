#!/bin/sh

# Comparing Metis and Scotch with varying number of cells and partitions
echo "Running benchmark 3"
echo "Running time is 2 minutes 13 seconds on a Sempron 3300+"
./benchmark metis ../../../data/meshes/torso.xml.gz 4 100000 metis_torso_4part.dat
./benchmark metis ../../../data/meshes/torso.xml.gz 16 100000 metis_torso_16part.dat
./benchmark metis ../../../data/meshes/torso.xml.gz 64 100000 metis_torso_64part.dat

./benchmark scotch ../../../data/meshes/torso.xml.gz 4 100000 scotch_torso_4part.dat
./benchmark scotch ../../../data/meshes/torso.xml.gz 16 100000 scotch_torso_16part.dat
./benchmark scotch ../../../data/meshes/torso.xml.gz 64 100000 scotch_torso_64part.dat
