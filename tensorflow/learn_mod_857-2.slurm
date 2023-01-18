#!/bin/bash -l
#SBATCH -J LEARN_2 
#SBATCH --nodes=1
#SBATCH --mem=118000
#SBATCH --constraint=HSW24
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=00:30:00

cd /scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
export OMP_NUM_THREADS=24
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc

export TF_XLA_FLAGS=--tf_xla_cpu_global_jit

outlog=/scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow

./learn2d.py 857-2 0 &> ${outlog}/learn_857-2.log 
# compute CNN 2D use parameters learned on 545-2
# srun --mpi=pmi2 -K1 -n ${numproc} --distribution=cyclic  ./SRoll4.py param545_ONE.py --onlyparam --rstep ${rstep} --cnn2d &> ${outlog}/test_TL_CNN2D.log


exit 0
