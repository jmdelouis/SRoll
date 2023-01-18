import numpy as np
import os

tab=[50,54,55]

os.system('mkdir -p RUN')

f = open('S857.sh','w')
for i in range(3):
    for itab in tab:
        os.system('sed "s;@XXX@;%03d;g;s;@XX@;%d;g" ../sroll/sroll3_857_data.troll_XXX.par > RUN/sroll3_857_sim.%03d_%d.par'%(i,itab,i,itab))
        os.system('sed "s;@XXX@;%03d;g;s;@XX@;%d;g" ../sroll/sroll3_857_XXX.qsub > RUN/sroll3_sim_%03d_%d.qsub'%(i,itab,i,itab))
        f.write("qsub RUN/sroll3_sim_%03d_%d.qsub\n"%(i,itab))

f.close()
