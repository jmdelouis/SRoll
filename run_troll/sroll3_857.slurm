#!/bin/bash -l
#SBATCH -J sroll3_857
#SBATCH --nodes=30
#SBATCH --mem=118000
#SBATCH --constraint=HSW24
#SBATCH --ntasks=256
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --time=10:30:00

cd /scratch/cnt0028/ias1717/SHARED/bware/sroll3_theo/extraball/SrollEx/sroll
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
export OMP_NUM_THREADS=1
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc
echo $SLURM_NTASKS
srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./troll_14tf /scratch/cnt0028/ias1717/SHARED/bware/sroll3_theo/extraball/SrollEx/run_troll/sroll3_857.troll.par &> /scratch/cnt0028/ias1717/SHARED/bware/sroll3_theo/extraball/SrollEx/run_troll/sroll3_857.log

exit 0
