## compute destriping
#mpirun -np 12 ./SRoll4.py param545_ONE.py --nocnn --rstep 10
#
## compute CNN 1D wloss=10^-2
#mpirun -np 12 ./SRoll4.py param545_ONE.py --isloss -2 --rstep 10
## compute CNN 1D wloss=10^-32
#mpirun -np 12 ./SRoll4.py param545_ONE.py --isloss -32 --rstep 10
## compute CNN 1D Learn weights on 545-2
#mpirun -np 12 ./SRoll4.py param545_ONE.py --doreference --rstep 10 --bolo '545-2'
## compute CNN 1D use parameters learned on 545-2
#mpirun -np 12 ./SRoll4.py param545_ONE.py --onlyparam --rstep 10

# compute CNN 2D wloss=10^-2
#mpirun -np 12 ./SRoll4.py param545_ONE.py --isloss -2 --rstep 10 --cnn2d --regul 10000
# compute CNN 2D wloss=10^-32
#mpirun -np 12 ./SRoll4.py param545_ONE.py --isloss -32 --rstep 10 --cnn2d --regul 10000
# compute CNN 2D Learn weights on 545-2
#mpirun -np 12 ./SRoll4.py param545_ONE.py --doreference --rstep 10 --cnn2d --bolo '545-1' --regul 10000
# compute CNN 2D use parameters learned on 545-2
mpirun -np 12 ./SRoll4.py param545_ONE.py --onlyparam --rstep 10 --cnn2d --regul 10000
mpirun -np 12 ./SRoll4.py param545_ONE.py --onlyparam --data --rstep 10 --cnn2d --regul 10000
#mpirun -np 12 ./SRoll4.py param545_ONE.py --data --rstep 10 --cnn2d --regul 10000 --isloss -32
