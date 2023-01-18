import numpy as np
from numpy import matrix
from numpy import linalg
import  matplotlib.pyplot as plt

f = open('matrix','r')
l = f.readlines()

data = []

for line in l:
  tmp_line = line.split("\t")
  tmp_line = tmp_line[0:len(tmp_line)-1]
  data.append(tmp_line)


m = matrix(data,dtype='double')  
print(m)
#print(m-m.transpose())

#m =linalg.inv(m)

plt.imshow(m,cmap ="jet")
plt.show()


 # ------------------------
vec = []
f1 = open('vec','r')
l1 = f1.readlines()

for line in l1:
  tmp_line = line.split("\t")
  for i in tmp_line[0:len(tmp_line)-1]:
    vec.append(i)


vec = np.array(vec,dtype='double')
print(vec.shape)

m = m+m[0,0]

r =linalg.solve(m,vec)
plt.plot(r)
plt.show()
