## compute destriping

./learn2d.py 857-1 -1 # save the difference once
./learn2d.py 857-2 0
./learn2d.py 857-3 0
./learn2d.py 857-4 0

mpirun -np 12 ./SRoll4.py param545_ONE.py --onlyparam --data --rstep 10 --regul 100 --bolo '857-1' --cnn2d --TFLEARN FSL_857-1_S-1_INITMODEL
mpirun -np 12 ./SRoll4.py param545_ONE.py --onlyparam --data --rstep 10 --regul 100 --bolo '857-2' --cnn2d --TFLEARN FSL_857-2_S0_INITMODEL
mpirun -np 12 ./SRoll4.py param545_ONE.py --onlyparam --data --rstep 10 --regul 100 --bolo '857-3' --cnn2d --TFLEARN FSL_857-3_S0_INITMODEL
mpirun -np 12 ./SRoll4.py param545_ONE.py --onlyparam --data --rstep 10 --regul 100 --bolo '857-4' --cnn2d --TFLEARN FSL_857-4_S0_INITMODEL

