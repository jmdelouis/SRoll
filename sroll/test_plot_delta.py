import numpy as np
from numpy import matrix
from numpy import linalg
import  matplotlib.pyplot as plt
import os
import sys

arg1 = int(sys.argv[1])
if(arg1 == 0):
  f = open('delta_old','r')
elif(arg1 == 1):
  f = open('delta','r')
elif(arg1 == 2):
  f = open('time','r')
elif(arg1 == 3):
  f = open('time_old','r')
elif(arg1 ==4):  
  f = open('time','r')
  f1 = open('time_old','r')

elif(arg1 ==5):
  f = open('delta','r')
  f1 = open('delta_old','r')

  

if(arg1==4 or arg1==5):
  l =f.readlines()
  l1 = f1.readlines()
  
  for i in l:
    data = i.split(";")
    data = data[0:-1]

    for i1 in l1:
      data1 = i1.split(";")
      data1 = data1[0:-1]

      data = np.array(data,dtype='double')
      data1 = np.array(data1,dtype='double')
      print(data)
      print(data1)
      plt.plot(data,color = 'red')
      plt.plot(data1,color = "blue")
      if(arg1==5):
        plt.yscale('log')
      plt.show()
      exit(0)
      

    
l = f.readlines()

data = []
color = ["red","blue","green","yellow","cyan","magenta","black","red","blue","green","yellow","cyan","magenta","black"]
it = 0

  
for line in l:
  line1 = line.split("//")
  for line in line1:
    if(l!="//"):
      tmp_line = line.split(";")
      print(tmp_line)
      data = tmp_line[1:len(tmp_line)-1]
      plop = np.array(data,dtype='double')
      plt.plot(plop,color = color[it])
      #plt.bar(np.arange(len(data)),plop)
      it= it+1

if(arg1<2): 
  plt.yscale('log')
plt.show()

