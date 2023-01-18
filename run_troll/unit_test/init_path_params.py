# script to init path for SRoll parametes file

import os 
import sys

   
#get paths    
print ("Path to SRoll:")
sroll_path = input( "> " )
print ("Path to data:")
data_path = input( "> " )
print ("Path for ouput:")
output_path = input( "> " )

# Read in the file
with open("template_sroll_unit_test_params.par", 'r') as file :
  filedata = file.read()

# Replace the target string
filedata = filedata.replace('/path_sroll/',str(sroll_path))
filedata = filedata.replace('/path_data/',str(data_path))
filedata = filedata.replace('/output_path/',str(output_path))

# Write the file out again
with open("sroll_unit_test_params.par", 'w') as file:
  file.write(filedata)


