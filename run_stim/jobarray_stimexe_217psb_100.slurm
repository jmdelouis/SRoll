#!/bin/bash -l

#SBATCH -J JAN18stim_217psb_100
#SBATCH --nodes=43
#SBATCH --constraint=HSW24
#SBATCH --ntasks=1024
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --mem=50GB
#SBATCH --time=3:0:00

#SBATCH --array=0-79%4

export MXM_TLS=self,shm,rc
cd /scratch/cnt0028/ias1717/SHARED/bware/sroll2r79/stim

PREFIX=/scratch/cnt0028/ias1717/SHARED/bware/sroll2r79/run_JAN18_217/JAN18

declare -a DETLIST=("217-5a" "217-5b" "217-6a" "217-6b" "217-7a" "217-7b" "217-8a" "217-8b")

# $(( )) does an arithmetic evaluation of its contents
DETECTOR="${DETLIST[$(($SLURM_ARRAY_TASK_ID % 8))]}"

# 8 detectors in DETLIST, 10 iterations per detector
FIRST_ITER=100
ITERNUM=$(($FIRST_ITER + ($SLURM_ARRAY_TASK_ID / 8) * 10))

# sample name: (...)/run_JAN18_217/JAN18_217-8b_190_217-8b.slurm.par
JOBNAME=${PREFIX}_${DETECTOR}_$(printf '%03d' $ITERNUM)_${DETECTOR}

srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./stimexe $JOBNAME.stim.par &> $JOBNAME.log

exit 0
