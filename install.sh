#!/bin/bash

# Create dirs
mkdir -p data

#Compilation of some files
g++ src/preprocess.cpp -O3 -march=native -fopenmp -std=c++11 -o src/preprocess.out
g++ util/get_denoising_stats.cpp -std=c++11 -O3 -march=native -o util/get_denoising_stats.out
g++ util/combine_into_fastq.cpp -std=c++11 -O3 -march=native -o util/combine_into_fastq.out
