#!/bin/sh
#PBS -q ice-s
#PBS -l select=1:ncpus=22:mpiprocs=22:mem=120gb
#PBS -N TEST_73
#PBS -o out_file
#PBS -j oe
cd /home/amitlal
mpiexec_mpt ./net7300 > net7300.log
