import numpy as np
import os

tab=[40,1240]

os.system('mkdir -p RUN')

f = open('D857.sh','w')
for itab in tab:
    os.system('sed "s;@XXX@;;g;s;@XX@;%d;g" ../sroll/sroll3_857_data.XXX.par > RUN/sroll3_data._%d.par'%(itab,itab))
    os.system('sed "s;@XXX@;;g;s;@XX@;%d;g" ../sroll/sroll3_857_data.qsub > RUN/sroll3_data_%d.qsub'%(itab,itab))
    f.write("qsub RUN/sroll3_data_%d.qsub\n"%(itab))

f.close()
