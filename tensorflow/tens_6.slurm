#!/bin/bash -l
#SBATCH -J sroll4
#SBATCH --nodes=6
#SBATCH --mem=118000
#SBATCH --constraint=HSW24
#SBATCH --ntasks=144
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=00:10:00

cd /scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
export OMP_NUM_THREADS=1
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc

export TF_XLA_FLAGS=--tf_xla_cpu_global_jit

srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./isroll2.py param545.py &> /scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow/testtens_6.log


exit 0
