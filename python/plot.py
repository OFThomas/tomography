#!/usr/local/bin/python3
#############################################################
# Title: Plot simulation data
#
# Date created: 19th June 2018
#
# Language:    Python 3
#
# Overview:    The program takes in a data file containing
#              simulation data and plots it. 
#
# Details:     The script performs the following steps:
#
# Usage: python3 linear-test.py
#
#############################################################

# Include
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import sys

if(len(sys.argv) != 2):
    print("Wrong number of command line args")
    print("You need to specify one data file to read")
    exit()

try:
    data_file = open(sys.argv[1], 'r')
except OSError as e:
    print(e)
    exit()

# Prepare variables to read in the data
#
av_distances = np.array([0,0,0])
non_physical = [] # Proportion of non-physical estimates
x_values = []

# Skip first 2 lines
for n in range(2) : next(data_file)

# Obtain simulation parameters
lineM = data_file.readline().split('=')
M = lineM[1].strip()
lineN = data_file.readline().split('=')
N = lineN[1].strip()
lineS = data_file.readline().split('=')
S = lineS[1].strip()

# Skip next 2 lines
for n in range(2) : next(data_file)

for line in data_file:
    current = line.split(',')
    if len(current) < 3 : break # at the end of the data block
    for k in range(len(current)):
        current[k] = current[k].strip()
    #av_distances[k][0:2] = current[0:2]
    #non_physical = current[3]
    x_values.append(current[0])
    av_distances = np.vstack((av_distances, current[1:4]))
    non_physical.append(current[4])

# Skip next 1 lines
#for n in range(1) : next(data_file)

# Find the total time line
while(True):
    lineT = data_file.readline().split('=')
    if len(lineT) == 2 : break
#print(lineT)
print(lineT)
total_time = lineT[1].strip()
#print("Time taken  =",total_time)

# Delete the first row of the array. This is a stupid workaround
# for the fact that the np array needs to start with a row in
# order to be able to add to it. There must be a better way than
# that
#
av_distances = np.delete(av_distances,0,0).astype(float)
non_physical = [float(i) for i in non_physical]
x_values = [float(i) for i in x_values]

print(av_distances)
print(non_physical)
print(x_values)

#print(av_distances[:,0])
#print(x_values)

# Plot the data
#
fig, ax1 = plt.subplots()
ax1.set_xlabel('Purity parameter')
ax1.set_ylabel('Estimate error distance')
ax1.plot(x_values, av_distances[:,0], '.',color='tab:red')
ax1.plot(x_values, av_distances[:,1], '.',color='tab:green')
ax1.plot(x_values, av_distances[:,2], '.',color='tab:blue')
ax1.legend(['Operator distance','Trace distance','Fidelity'])

ax2=ax1.twinx()
ax2.set_ylabel('Probability of non-physical estimate')
ax2.plot(x_values,non_physical, '+', color='tab:brown')
atext = AnchoredText("Number of purity parameters: " + str(M) + "\n"
                     +"Density matrices: "+ str(N) + "\n"
                     +"Samples per measurement: " + str(S) + "\n"
                     +"Total running time:  " + str(total_time),
                     loc=2)
ax2.add_artist(atext)

fig.tight_layout()
plt.show()
