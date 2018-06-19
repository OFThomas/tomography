#!/usr/local/bin/python3
#############################################################
# Title: Testing the linear estimator
#
# Date created: 19th June 2018
#
# Language:    Python 3
#
# Overview:    The point of this program is to test the
#              reliability of the linear estimator for
#              density matrices with different levels of
#              purity. 
#
# Details:     The script performs the following steps:
#
#              1) Set x in [0,1] as the purity parameter;
#                 pick a random unitary U; and compute
#
#
#                          -       -
#                         | x     0 |
#                 p = U * |         | * U^
#                         | 0   1-x |
#                          -       -
#
#                 The closer x is to 1 or 0, the more
#                 pure is the state.
#
#
#              2) Fix a set of measurement operators,
#                 and generate sample data from the
#                 distribution of outcomes implied by
#                 the state p.
#
#              3) Compute the density matrix using the
#                 linear estimator. Then compute the
#                 distance between the estimate and the
#                 true density matrix for each of the
#                 different distances (operator norm,
#                 Hilbert-Schmidt norm, fidelity).
#
#              4) Repeat steps 1-3 for different values
#                 of x and repeat each value of x a large
#                 number of times. Store all the distances.
#
#              5) Perform averages of each distance at each
#                 value of x. Plot the average distances
#                 as a function of x for all the values of
#                 x.
#
# Usage: python3 linear-test.py
#
#############################################################

# Include
import importlib
import numpy as np
from scipy.stats import unitary_group as ug
import scipy as sc

################
## SIMULATION ##
################

# Step 1: Prepare the density matrix
#
# The purity parameter x is picked between 0
# and 1.
#
# Note: any time a numerical check is performed
# and printed out, I've rounded the result to
# make it more readable. Set the decimal places
# to keep using the dp variable.
#
dp = 5
x = 0.6
import simulation
importlib.reload(simulation)
dens = simulation.random_density(x,'none',dp)


