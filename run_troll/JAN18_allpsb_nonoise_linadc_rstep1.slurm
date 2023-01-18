#!/bin/bash -l
#SBATCH -J JAN18_allpsb_nonoise_linadc_rstep1
#SBATCH --nodes=43
#SBATCH --constraint=HSW24
#SBATCH --ntasks=1024
#SBATCH --ntasks-per-node=24
#SBATCH --threads-per-core=1
#SBATCH --mem=50GB
#SBATCH --time=8:0:00

cd /scratch/cnt0028/ias1717/SHARED/bware/sroll2r94/sroll
# source ./srollex_setenv.sh is not needed because environment variable are propagated to the job
# use Mellanox reliable connection http://www.mellanox.com/related-docs/prod_software/Mellanox_MXM_README_v2.0.pdf
export MXM_TLS=self,shm,rc
srun --mpi=pmi2 -K1 -n $SLURM_NTASKS --distribution=cyclic ./troll /scratch/cnt0028/ias1717/SHARED/bware/sroll2r94/run_troll/JAN18_allpsb_nonoise_linadc_rstep1.troll.par &> /scratch/cnt0028/ias1717/SHARED/bware/sroll2r94/run_troll/JAN18_allpsb_nonoise_linadc_rstep1.log

# enable the following line to automatically merge fits IQU after sroll. IT CONSUMES CPU TIME
#srun --mpi=pmi2 -K1 -n 43 --distribution=cyclic /scratch/cnt0028/ias1717/SHARED/bware/sroll2r94/run_troll/merge_fits_maps.py /scratch/cnt0028/ias1717/SHARED/bware/sroll3test_MAP/JAN18_allpsb_nonoise_linadc_rstep1*I.fits >>/scratch/cnt0028/ias1717/SHARED/bware/sroll2r94/run_troll/JAN18_allpsb_nonoise_linadc_rstep1.log 2>&1

exit 0
