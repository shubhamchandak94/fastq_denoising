#!/bin/bash

# Create dirs
mkdir -p data

#Compilation of some files
g++ src/preprocess.cpp -O3 -march=native -fopenmp -std=c++11 -o src/preprocess.out
