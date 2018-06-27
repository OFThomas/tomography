#!/usr/local/bin/python3
#############################################################
# Title: Quantum state tomography simulation
#
# Date created: 6th June 2018
#
# Language:    Python 3
#
# Overview:    Simulation of linear quantum state
#              tomography on randomly generated
#              1 qubit states. Linear quantum
#              state tomography involves estimating
#              the density matrix of the 1 qubit
#              state by measuring successive
#              instances of the state in an
#              informationally complete basis of
#              operators, and then using the relative
#              frequencies of the results to determine
#              the coefficients of the density matrix
#              expressed in this basis of operators:
#
#              E.g. (p = density matrix)
#
#              p = tr(pI)I + tr(pX)X + tr(pY)Y + tr(pZ)Z
#
#              Since tr(pS) is the expectation value
#              of S in the state p, repeated measurement
#              of S followed by averaging will result in
#              tr(pS). Applying the proceedure to S = X,
#              Y, and Z will yield three of the four
#              coefficients. The fourth is obtained by
#              the condition that tr(p) = 1.
#
# Details:     The script performs the following steps:
#
#              1) Generate random one qubit density
#                 matrices as follows: pick x in [0,1];
#                 pick a random unitary U; and compute
#
#
#                          -       -
#                         | x     0 |
#                 p = U * |         | * U^
#                         | 0   1-x |
#                          -       -
#
#                 x should be picked in such a way
#                 that the resulting distribution for
#                 p should be suitably uniform (caveat!)
#                 Alternatively, it should be possible
#                 for the user to put in a specific
#                 density matrix.
#
#              2) Fix a set of measurement operators,
#                 and generate sample data from the
#                 distribution of outcomes implied by
#                 the state p.
#
#              3) Choose an estimator for p, and
#                 compute the estimated state p~
#                 using the data generated in step
#                 2. Then compare p~ with p.
#
#              4) Compare the generated density matrix
#                 with the original density matrix
#                 using some measure of distance (e.g.
#                 Hilbert-Schmidt norm, operator norm,
#                 etc.)
#
# Notes:       1) Throughout, the Hermitian conjugate
#                 of a matrix A is is denoted A^
#
# Usage: python3 state-tomography.py
#
#############################################################

# Include
import importlib
import numpy as np
from scipy.stats import unitary_group as ug
import scipy as sc
from input import get_user_value

################
## SIMULATION ##
################

# Step 1: Prepare the density matrix
#
# The variable x is drawn from a uniform
# distribution on [0,1]. The random
# unitary group U is is generated using
# scipy.stats.unitary_group.rvs
#
# Note: any time a numerical check is performed
# and printed out, I've rounded the result to
# make it more readable. Set the decimal places
# to keep using the dp variable.
#
dp = 5
import simulation
importlib.reload(simulation)
x = np.random.uniform(0,1) # Generate x
print("Picked eigenvalues:\n\n\t ",x,",",1-x)
dens = simulation.random_density(x)
print("\nThe random density matrix is\n\n")
print(dens)
    
# Step 2: Generate measurement data
#
# This step uses the following single
# qubit measurements: X, Y and Z. Each
# measurement results in either +1 or
# -1, with a distribution determined
# by the density matrix.
#
# Suppose X is measured. The eigenvectors
# of X are |+1> and |-1>, corresponding
# to the outcomes of the measurement
# (+1 or -1). The probability of each
# possibility in the state p is given by
# tr(p|+1><+1|) or tr(p|-1><-1|),
# depending on whether the probability
# of +1 or -1 is required. If |+1> is
# expressed as a vector [a,b], then
# <+1| is expressed as [a,b]^T, and
# |+1><+1| = [a,b][a,b]^T. 
#
I = np.matrix([[1,0],[0,1]])
X = np.matrix([[0,1],[1,0]])
Y = np.matrix([[0,-1j],[1j,0]])
Z = np.matrix([[1,0],[0,-1]])
measurements = np.array([X,Y,Z,I])
meas_ops = {'X':X, 'Y':Y, 'Z':Z}
print("\n===== Generate measurement data =====\n")
print("The X measurements are simulated by")
print("computing the projectors P = |+><+| and ")
print("P = |-><-| for each measurement outcome")
print("and then computing the probability using")
print("prob = tr(pP), where p is the density")
print("matrix. Then samples are taken from")
print("a discrete distribution defined by the")
print("measurement probabilities. The process")
print("is repeated for Y and Z.")

# Preliminaries: compute the projectors

I = np.matrix([[1,0],[0,1]])
X = np.matrix([[0,1],[1,0]])
Y = np.matrix([[0,-1j],[1j,0]])
Z = np.matrix([[1,0],[0,-1]])

values_X, vectors_X = np.linalg.eig(X)
proj_X = np.zeros([2,2,2])
proj_X[0,:,:] = np.matmul(vectors_X[:,0], np.matrix.getH(vectors_X[:,0]))
proj_X[1,:,:] = np.matmul(vectors_X[:,1], np.matrix.getH(vectors_X[:,1]))
print("\nThe projector for",values_X[0],"outcome of X is\n\n",proj_X[0,:,:])
print("\nThe projector for",values_X[1],"outcome of X is\n\n",proj_X[1,:,:])


values_Y, vectors_Y = np.linalg.eig(Y)
proj_Y = np.zeros([2,2,2],dtype='complex')
proj_Y[0,:,:] = np.matmul(vectors_Y[:,0], np.matrix.getH(vectors_Y[:,0]))
proj_Y[1,:,:] = np.matmul(vectors_Y[:,1], np.matrix.getH(vectors_Y[:,1]))
print("\nThe projector for",values_Y[0],"outcome of Y is\n\n",proj_Y[0,:,:])
print("\nThe projector for",values_Y[1],"outcome of Y is\n\n",proj_Y[1,:,:])


values_Z, vectors_Z = np.linalg.eig(Z)
proj_Z = np.zeros([2,2,2])
proj_Z[0,:,:] = np.matmul(vectors_Z[:,0], np.matrix.getH(vectors_Z[:,0]))
proj_Z[1,:,:] = np.matmul(vectors_Z[:,1], np.matrix.getH(vectors_Z[:,1]))
print("\nThe projector for",values_Z[0],"outcome of Z is\n\n",proj_Z[0,:,:])
print("\nThe projector for",values_Z[1],"outcome of Z is\n\n",proj_Z[1,:,:])

S = get_user_value("\nChoose the number of measurements in each basis","integer")  

X_data = simulation.simulate(dens,proj_X,values_X,S)
Y_data = simulation.simulate(dens,proj_Y,values_Y,S)
Z_data = simulation.simulate(dens,proj_Z,values_Z,S)

################
## ESTIMATION ##
################

# Step 3: Estimate p using the simulated data
#
# The estimator here involves computing the
# mean of data sets, and using them as
# estimates for tr(pX), tr(pY) and tr(pZ)
#
# Then tr(pI) is computed by requiring that
# the density matrix be normalised
#
import estimation
importlib.reload(estimation)
dens_est = estimation.linear_estimate_XYZ(X_data, Y_data, Z_data)

print("The estimate for p is:\n\n",dens_est,"\n")
print("The original density matrix was:\n\n", dens,"\n")

# Step 4: Compute the distance between p and p~
#
# Compute the distance between dens and
# dens_est. Distance could be computed
# using the metric that arises from the
# operator norm, Hilbert-Schmidt norm
# also called the Frobenius norm), etc.
#
# The operator norm of A is
#
#     ||A|| = supp |Ax|
#            |x|=1
#
# The Hilbert-Schmidt norm of A is:
#
#     ||A|| = sqrt[ tr(A^ A) ]
#
# A metric d is always obtained from
# a norm using
#
#     d(A,B) = ||A - B||
#
# The fidelity distance is defined as
# follows
#
#     F(A,B) = tr[ sqrt(sqrt(A) B sqrt(A)) ]
#     d(A.B) = arccos[F(A,B)]
#
# The np.norm function implements
# several norms depending on the
# value of the second parameter:
#
#     'fro' = Frobenius norm
#
from stats import * 
print("======================= Summary statistics =======================\n")
print("The number of simulated samples for each measurement was\u001b[36m",X_data.size,"\u001b[37m\n")
variances = {}
# Get the purity of the density matrix estimate
eigenvalues, eigenvectors = np.linalg.eig(dens_est)
# Assume eigenvalues are real
print("The eigenvalues of the estimated density matrix are:\n\n\t",
      np.around(eigenvalues.real[0],dp), "and",
      np.around(eigenvalues.real[1],dp))
print("\nThe purity parameter of the estimated density matrix is",
      np.maximum(np.around(eigenvalues.real[0],dp),
      np.around(eigenvalues.real[1],dp)))
print("\nThe variances in the simulated data are:\n")

print("\tThe variance in the X samples was",np.var(X_data))
print("\tThe variance in the Y samples was",np.var(Y_data))
print("\tThe variance in the Z samples was",np.var(Z_data))

print("\nDistances between the original density matrix and the estimate:")
print("\n\tIn the operator norm:\t\t", distance_op(dens, dens_est))
print("\tIn the Hilbert Schmidt norm:\t",distance_trace(dens, dens_est))
print("\tFidelity distance:\t\t",distance_fid(dens, dens_est))
print()
