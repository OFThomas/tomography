#!/usr/local/bin/python3
#############################################################
# Title: Linear adaptive test
#
# Date created: 19th June 2018
#
# Language:    Python 3
#
# Overview:    The program uses a simple adaptive strategy
#              to attempt to improve the linear estimate.
#              The first part of the data is used to choose
#              a measurement basis, which is then used for
#              all subsequent measurements.
#
# Details:     The script performs the following steps:
#
# Usage: python3 linear-adaptive.py
#
#############################################################

# Include
import importlib
import numpy as np
from numpy import sin, cos
from scipy.stats import unitary_group as ug
import scipy as sc
import simulation
importlib.reload(simulation)
import estimation
importlib.reload(estimation)
import stats
importlib.reload(stats)
import estimation
importlib.reload(estimation)
import cProfile
import pstats
from progress import *

pr = cProfile.Profile()
pr.enable()

# ======= Test parameter ===============================
M = 2000  # Number of purity parameters x to try
x_start = 0 # Specify purity parameter x range
x_end = 1
N = 500  # Number of random density matrices per x value
S = 500  # Number of samples of each measurement to
         # simulate for each density matrix
S_1 = 50 # Number of samples of X, Y and Z to use to
         # perform the prelimiary estimate to adapt
         # the basis.
# ======================================================

# Seed random number generator
np.random.seed()

av_distances = np.zeros([M,3])
non_physical = np.zeros(M) # Proportion of non-physical estimates

# Preliminaries: compute the projectors

I = np.matrix([[1,0],[0,1]])
X = np.matrix([[0,1],[1,0]])
Y = np.matrix([[0,-1j],[1j,0]])
Z = np.matrix([[1,0],[0,-1]])

# This is a terrible function
proj_X, proj_Y, proj_Z, values_X, values_Y, values_Z = simulation.projectors(X,Y,Z)

# Open a file for writing
file = open("linear_adaptive_1_python.dat", "w")
file.write("Distances between estimated and original density matrices using various distances:\n\n")
file.write("Number of purity values tried = "+str(M)+"\n")
file.write("Number of density matrices per purity parameter = "+str(N)+"\n")
file.write("Total number of measurements for each of X, Y and Z = "+str(S)+"\n\n");
file.write("PURITY, \tOPERATOR, \tTRACE, \t\tFIDELITY, \tNON PHYSICAL\n")

# Define x -- put a loop here ----------------------------------------- LOOP for x between 0 and 1
#
# This loop runs through different values of the purity parameter x,
# and tests the ability of the linear estimator in each case
#

dist_op = np.zeros([N,1])
dist_trace = np.zeros([N,1])
dist_fid = np.zeros([N,1])

dp = 5

# Rotate v about unit vector u
# by angle p
def axis_rot(u,v,p) :
    # Rotation of one vector about unit axis (u_x, u_y, u_z) by angle p
    r11 = cos(p) + u[0]**2 * (1 - cos(p))
    r22 = cos(p) + u[1]**2 * (1 - cos(p))
    r33 = cos(p) + u[2]**2 * (1 - cos(p))
    r12 = u[0] * u[1] * (1 - cos(p)) - u[2] * sin(p)
    r21 = u[0] * u[1] * (1 - cos(p)) + u[2] * sin(p)
    r13 = u[0] * u[2] * (1 - cos(p)) - u[1] * sin(p)
    r31 = u[0] * u[2] * (1 - cos(p)) + u[1] * sin(p)
    r23 = u[1] * u[2] * (1 - cos(p)) - u[0] * sin(p)
    r32 = u[1] * u[2] * (1 - cos(p)) + u[0] * sin(p)
    
    Ru = np.array([[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]])

    # Get the other two measurement directions
    w = np.matmul(Ru, v)

    return w

#for k,n in itertools.product(range(M),range(N)):
for k in range(M):
    non_physical_count = 0 # Temporary counter for non-physical estimates
    
    # Loop N times for each value of x ------------------ inner loop -- N trials for each x
    #
    # This loop generates N random density matrices for each fixed value of x
    # which used to simulate measurement data and run the estimator
    #
    for n in range(N):
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
        x = x_start + k * (x_end - x_start)/M
        #pr.enable()
        dens = simulation.random_density(x)
        #values_dens,vectors_dens = np.linalg.eig(dens)
        #pr.disable()

        # Step 2: Generate measurement data
        #
        # Generate data for X, Y and Z measurements.
        # Only generate S_1 samples for each measurement.
        # These will be used to estimate a preliminary
        # density matrix, that will be used to adapt the
        # basis.
        #
        X_data = simulation.simulate(dens,proj_X,values_X,S_1)
        Y_data = simulation.simulate(dens,proj_Y,values_Y,S_1)
        Z_data = simulation.simulate(dens,proj_Z,values_Z,S_1)

        # Step 3: Estimate the preliminary density matrix
        #
        # Compute linear estimator
        #
        # Then tr(pI) is computed by requiring that
        # the density matrix be normalised
        #
        dens_est = estimation.linear_estimate_XYZ(X_data, Y_data, Z_data)

        # Step 4: Adapt the measurement basis
        #
        # First, find the Bloch vector for the density matrix.
        # This is obtained using the Hilbert-Schmidt inner
        # product (A,B) = Tr(A^ B):
        #
        # V = 1/2 * [Tr(dens_est * X),
        #            Tr(dens_est * Z),
        #            Tr(dens_est * Z)]
        #
        # Then, rotate the vector V in an arbitrary direction by
        # 0.95531 rad (54 degrees).
        #
        # The new Bloch vector W will be one of the directions for the
        # new measurements. Then rotate this new vector around the first
        # Bloch vector by 120 degrees and then 240 degrees. These vectors
        # will be the other two measurement Bloch vectors.
        #
        # Use the Bloch vectors to obtain Hermitian matrices using the
        # expansion H = (1/2)(I + V.X) where V is the Bloch vector and
        # X is the vector of Pauli matrices.
        #
        V = np.array([np.trace(np.matmul(dens_est,X)),
                      np.trace(np.matmul(dens_est,Y)),
                      np.trace(np.matmul(dens_est,Z))]).transpose()
        theta = np.pi /2#np.arccos(1/np.sqrt(3))
        R1 = np.array([[1,0,0],[0,cos(theta),-sin(theta)],[0,sin(theta),cos(theta)]])
        W1 = np.matmul(R1,V)
        p = (2 * np.pi)/3

        # Axis of rotation
        u = V / np.linalg.norm(V)
        print(np.linalg.norm(u))
        
        # Obtain the other measurement axes
        W2 = axis_rot(u, W1, p)
        W3 = axis_rot(u, W1, 2*p)
        W4 = axis_rot(u, W1, 3*p)
        print(W1)
        print(W2)
        print(W3)
        print(W4)
        
        # Check inner products between basis Bloch vectors
        # They should be orthonormal
        print("Here", np.dot(W1,W2))
        print("Here", np.dot(W1,W3))
        print("Here", np.dot(W2,W3))
        exit()
        # Generate the measurement matrices
        M1 = (1/2) * (I + W1[0]*X + W1[1]*Y + W1[2]*Z)
        M2 = (1/2) * (I + W2[0]*X + W2[1]*Y + W2[2]*Z)
        M3 = (1/2) * (I + W3[0]*X + W3[1]*Y + W3[2]*Z)

        # Generate projectors
        proj_M1, proj_M2, proj_M3 = simulation.projectors(M1,M2,M3)

        # Step 5: Simulate new measurements in the new basis
        #
        # This step simulates data in the new measurement
        # basis. These new measurements are then used to
        # obtain the final estimate of the density matrix.
        #
        M1_data = simulation.simulate(dens,proj_X,values_X,S_1)
        M2_data = simulation.simulate(dens,proj_Y,values_Y,S_1)
        M3_data = simulation.simulate(dens,proj_Z,values_Z,S_1)
        
        
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
        eigenvalues = np.linalg.eigvals(dens_est)
        if eigenvalues[0] < 0 or eigenvalues[1] < 0:
            non_physical_count = non_physical_count + 1
        
    # Step 5: Average the distances 
    #
    # Average the distances for each value of x
    #
    av_distances[k,:] = [np.mean(dist_op), np.mean(dist_trace), np.mean(dist_fid)]
    non_physical[k] = non_physical_count/N
    p = (k+1)/M
    show_progress(pr,p)
    file.write("{0:.5f},\t{1:.5f},\t{2:.5f}, \t{3:.5f}, \t{4:.5f}\n".format(x, np.mean(dist_op),
                                                                            np.mean(dist_trace),
                                                                            np.mean(dist_fid),
                                                                            non_physical[k]))


pr.disable()
ps = pstats.Stats(pr)
total_time = ps.total_tt

file.write("\nTotal running time = "+str(np.around(total_time,3))+"s\n")    
file.close
