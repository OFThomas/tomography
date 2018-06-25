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
#                 linear estimator.
#
#              4) Compute the distance between the
#                 estimate and the true density matrix
#                 for each of the different distances
#                 (operator norm, Hilbert-Schmidt norm,
#                 fidelity).
#
#              5) Perform averages of each distance at each
#                 value of x. Plot the average distances
#                 as a function of x for all the values of
#                 x.
#
#              6) Repeat steps 1-4 for different values
#                 of x and repeat each value of x a large
#                 number of times. Store all the distances.
#
# Usage: python3 linear-test.py
#
#############################################################

# Include
import importlib
import numpy as np
from scipy.stats import unitary_group as ug
import scipy as sc
import simulation
importlib.reload(simulation)
import estimation
importlib.reload(estimation)
import stats
importlib.reload(stats)
import matplotlib.pyplot as plt
import estimation
importlib.reload(estimation)
from matplotlib.offsetbox import AnchoredText
import cProfile
import pstats

pr = cProfile.Profile()
pr.enable()

# ======= Test parameter ===============================
M = 20 # Number of purity parameters x to try
x_start = 0 # Specify purity parameter x range
x_end = 1
N = 50   # Number of random density matrices per x value
S = 50   # Number of samples of each measurement to
         # simulate for each density matrix 
# ======================================================

av_distances = np.zeros([M,3])
non_physical = np.zeros([M,1]) # Proportion of non-physical estimates

# Define x -- put a loop here ----------------------------------------- LOOP for x between 0 and 1
#
# This loop runs through different values of the purity parameter x,
# and tests the ability of the linear estimator in each case
#

for k in range(0,M):
    dist_op = np.zeros([N,1])
    dist_trace = np.zeros([N,1])
    dist_fid = np.zeros([N,1])
    non_physical_count = 0 # Temporary counter for non-physical estimates
    
    # Loop N times for each value of x ------------------ inner loop -- N trials for each x
    #
    # This loop generates N random density matrices for each fixed value of x
    # which used to simulate measurement data and run the estimator
    #
    for n in range(0,N):
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
        x = x_start + k * (x_end - x_start)/M
        dens = simulation.random_density(x)
        
        # Step 2: Generate measurement data
        #
        # Generate data for X, Y and Z measurements
        #
        I = np.matrix([[1,0],[0,1]])
        X = np.matrix([[0,1],[1,0]])
        Y = np.matrix([[0,-1j],[1j,0]])
        Z = np.matrix([[1,0],[0,-1]])
        measurements = np.array([X,Y,Z,I])
        meas_ops = {'X':X, 'Y':Y, 'Z':Z}
        sim_dat = simulation.simulate(dens,meas_ops,S)
        
        # Step 3: Estimate density matrix
        #
        # Compute linear estimator
        #
        # Then tr(pI) is computed by requiring that
        # the density matrix be normalised
        #
        dens_est = estimation.linear_estimate_XYZ(sim_dat)
        
        # Step 4: Compute and the distances
        #
        # Compute distances between the estimated
        # and true density matrix using the
        # different distance fuctions.
        #
        dist_op[n] = stats.distance_op(dens, dens_est)
        dist_trace[n] = stats.distance_trace(dens, dens_est)
        dist_fid[n] = stats.distance_fid(dens, dens_est)

        # Count the number of non-physical matrices
        #
        eigenvalues, eigenvectors = np.linalg.eig(dens_est)
        if eigenvalues[0] < 0 or eigenvalues[1] < 0:
            non_physical_count = non_physical_count + 1 
    # Step 5: Average the distances 
    #
    # Average the distances for each value of x
    #
    print(eigenvalues[0])
    av_distances[k,:] = [np.mean(dist_op), np.mean(dist_trace), np.mean(dist_fid)]
    non_physical[k] = non_physical_count/N

pr.disable()
pr.create_stats()

ps = pstats.Stats(pr)
total_time = ps.total_tt
pr.print_stats()


# Step 6: Process the data 
#
# Store in a file, plot, etc.
x_values = np.linspace(x_start, x_end, M)

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
                     +"Total running time:  " + str(np.around(total_time,3)) + "s",
                     loc=2)
ax2.add_artist(atext)

fig.tight_layout()
plt.show()
