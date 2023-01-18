#!/bin/bash -l
#SBATCH -J sroll4
#SBATCH --nodes=24
#SBATCH --mem=118000
#SBATCH --constraint=HSW24
#SBATCH --ntasks=576
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=01:00:00

cd /scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
export OMP_NUM_THREADS=1
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc

export TF_XLA_FLAGS=--tf_xla_cpu_global_jit

srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./isroll.py param545.py &> /scratch/cnt0028/ias1717/SHARED/bware/sroll22_trick/SrollEx/tensorflow/testtens_24.log

# enable the following line to automatically merge fits IQU after sroll. IT CONSUMES CPU TIME
#srun --mpi=pmi2 -K1 -n 43 --distribution=cyclic /scratch/cnt0028/ias1717/SHARED/bware/wrk_sroll21/run_stim/merge_fits_maps.py /scratch/cnt0028/ias1717/smottet/s3_sroll21_MAP/s21sim3_000_353ghz*I.fits >>/scratch/cnt0028/ias1717/SHARED/bware/wrk_sroll21/run_stim/s21sim3_000_353ghz.log 2>&1

exit 0
