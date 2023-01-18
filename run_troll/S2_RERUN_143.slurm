#!/bin/bash -l
#SBATCH -J S2_RERUN_143
#SBATCH --constraint=HSW24
#SBATCH --nodes=86
#SBATCH --ntasks=2048
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --mem=50GB
#SBATCH --time=02:00:00

cd /scratch/cnt0028/ias1717/SHARED/bware/sroll2r99_jmd/extraball/SrollEx/sroll
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc
srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./troll /scratch/cnt0028/ias1717/SHARED/bware/sroll2r99_jmd/extraball/SrollEx/run_troll/param143_data_sroll2r79_dipHFI17.txt &> /scratch/cnt0028/ias1717/SHARED/bware/sroll2r99_jmd/extraball/SrollEx/run_troll/S2_RERUN_143.log

exit 0
