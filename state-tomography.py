#############################################################
# Title: Quantum state tomography simulation
#
# Language:    Python
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
#                 p = U * |         | * U^+
#                         | 0   1-x |
#                          -       -
#                            -      -
#                           | a       |
#                 where U = |        |
#                           |        |
#                            -      -
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
#                 2. Compare p~ with p.
#
# Date created: 6th June 2018
#
#############################################################

import numpy as np
from scipy.stats import unitary_group as ug

# Step 1: Prepare the density matrix
#
# Note: any time a numerical check is performed
# and printed out, I've rounded the result to
# make it more readable. Set the decimal places
# to keep using the dp variable.
#
dp = 5
x = np.random.uniform(0,1) # Generate x
print("The value of x is: ", x, "\n")
realMat = np.random.random((2,2))
U = ug.rvs(2) # Random unitary 
print("The random unitary is: \n")
print(U,"\n")
# Check that U is actually unitary
U_dag = np.matrix.getH(U)
test = np.matmul(U_dag,U)
print("Check that U * U^+ = I: \n")
print(np.around(test,dp),"\n")
# Compute the density matrix
diag = np.matrix([[x,0],[0,1-x]])
dens = np.matmul(U,np.matmul(diag,U_dag))
print("The density matrix is:\n\n", dens,"\n")
# Check the density matrix
print("The trace of the density matrix is:", np.around(np.matrix.trace(dens),dp))
print("The eigenvalues are:", np.around(np.linalg.eig(dens)[0],dp), "which should both be positive.")

# Step 2: Generate measurement data
#
# 
#
#
#
#
