#!/bin/bash
#PBS -q mpi
#PBS -l select=1:ncpus=28:mem=60g
#PBS -l walltime=01:00:00
#PBS -m n

cd /home1/datahome/tfoulqui/MySroll/extraball/SrollEx/
source srollex_setenv.sh
which python

cd sroll

mpirun -n 8 ./troll_14tf /home1/datahome/tfoulqui/MySroll/extraball/SrollEx/run_troll/sroll3_857_data.troll_jmd.par


