#!/bin/bash -l
#SBATCH -J 857_2
#SBATCH --nodes=10
#SBATCH --mem=118000
#SBATCH --constraint=HSW24
#SBATCH --ntasks=240
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=00:30:00

cd /scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
export OMP_NUM_THREADS=1
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc

export TF_XLA_FLAGS=--tf_xla_cpu_global_jit

numproc=240
rstep=1
outlog=/scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow/

# compute destriping
srun --mpi=pmi2 -K1 -n ${numproc} --distribution=cyclic ./SRoll4.py param545_ONE.py --data --rstep ${rstep} --bolo '857-2' --nocnn &> ${outlog}/ref_fsl_857-2.log


exit 0
