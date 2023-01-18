#!/bin/bash
#PBS -q mpi_4
#PBS -l select=4:ncpus=28:mem=60g
#PBS -l walltime=06:30:00
#PBS -N Sroll_NORM

source /usr/share/Modules/3.2.10/init/bash

module purge
module load impi/5.1.3.258
export OMP_NUM_THREADS=1

cd ~/sroll4/extraball/SrollEx/tensorflow

numproc=64
rstep=100

# compute destriping
mpirun -np ${numproc} ./SRoll4.py param545_ONE.py --nocnn --rstep ${rstep} &> /scratch/home1/jmdeloui/test_TL_NORM.log

# compute CNN 1D wloss=10^-2
mpirun -np ${numproc}  ./SRoll4.py param545_ONE.py --isloss -2 --rstep ${rstep} &> /scratch/home1/jmdeloui/test_TL_CNN1D_2.log
# compute CNN 1D wloss=10^-32
mpirun -np ${numproc} ./SRoll4.py param545_ONE.py --isloss -32 --rstep ${rstep} &> /scratch/home1/jmdeloui/test_TL_CNN1D_32.log
# compute CNN 1D Learn weights on 545-2
mpirun -np ${numproc} ./SRoll4.py param545_ONE.py --doreference --rstep ${rstep} --bolo '545-2' &> /scratch/home1/jmdeloui/test_TL_CNN1D_REF.log
# compute CNN 1D use parameters learned on 545-2
mpirun -np ${numproc} ./SRoll4.py param545_ONE.py --onlyparam --rstep ${rstep} &> /scratch/home1/jmdeloui/test_TL_CNN1D.log

# compute CNN 2D wloss=10^-2
mpirun -np ${numproc} ./SRoll4.py param545_ONE.py --isloss -2 --rstep ${rstep} --cnn2d &> /scratch/home1/jmdeloui/test_TL_CNN2D_2.log
# compute CNN 2D wloss=10^-32
mpirun -np ${numproc} ./SRoll4.py param545_ONE.py --isloss -32 --rstep ${rstep} --cnn2d &> /scratch/home1/jmdeloui/test_TL_CNN2D_32.log
# compute CNN 2D Learn weights on 545-2
mpirun -np ${numproc}  ./SRoll4.py param545_ONE.py --doreference --rstep ${rstep} --cnn2d --bolo '545-2' &> /scratch/home1/jmdeloui/test_TL_CNN2D_REF.log
# compute CNN 2D use parameters learned on 545-2
mpirun -np ${numproc}  ./SRoll4.py param545_ONE.py --onlyparam --rstep ${rstep} --cnn2d &> /scratch/home1/jmdeloui/test_TL_CNN2D.log



