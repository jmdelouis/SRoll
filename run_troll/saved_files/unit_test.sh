
#UNIT TEST -- Sroll

#!/bin/bash
#PBS -q mpi
#PBS -l select=1:ncpus=28:mem=60g
#PBS -l walltime=01:00:00
#PBS -m n

HOST=$( hostname )

case "${HOST}" in
  login[0-9]*)
    echo "Run on OCCIGEN"
    ;;
  datarmor[0-9]*)
    echo "Run on DATARMOR"
    cd /home1/datahome/tfoulqui/MySroll/extraball/SrollEx/
    source srollex_setenv.sh
    cd sroll
    export OMP_NUM_THREADS=1
    mpirun -n 32 ./troll_14tf /home1/datahome/tfoulqui/MySroll/extraball/SrollEx/run_troll/sroll3_857.troll.par --rstep 100
    ;;
  *)
    echo "WARNING: unknown host '${HOST}'! (Nothing done)"
    return 1
    ;;
esac
