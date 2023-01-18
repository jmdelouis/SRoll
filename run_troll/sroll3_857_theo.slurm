#!/bin/bash
#PBS -q mpi
#PBS -l select=1:ncpus=28:mem=60g
#PBS -l walltime=01:00:00
#PBS -m n

cd /home1/datahome/tfoulqui/MySroll/extraball/SrollEx/sroll
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
export OMP_NUM_THREADS=1
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc
srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./troll_14tf /home1/datahome/tfoulqui/MySroll/extraball/SrollEx/run_troll/sroll3_857_theo.troll.par &> /home1/datahome/tfoulqui/MySroll/extraball/SrollEx/run_troll/sroll3_857_theo.log

exit 0
