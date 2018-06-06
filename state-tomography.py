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
# Date created: 6th June 2018
#
#############################################################

# Include
import numpy as np
from scipy.stats import unitary_group as ug

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
dens = U * diag * U_dag
print("The density matrix is:\n\n", dens,"\n")
# Check the density matrix
print("The trace of the density matrix is:", np.around(np.matrix.trace(dens),dp))
print("The eigenvalues are:", np.around(np.linalg.eig(dens)[0],dp), "which should both be positive.")

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
X = np.matrix([[0,1],[1,0]])
Y = np.matrix([[0,-1j],[1j,0]])
Z = np.matrix([[1,0],[0,-1]])
print("The measurements are X, Y and Z:\n\n",X,"\n\n",Y,"\n\n and \n\n",Z,".\n")
# Compute the eigenvectors and eigenvalues of X, Y and Z
X_values, X_vectors = np.linalg.eig(X)
Y_values, Y_vectors = np.linalg.eig(Y)
Z_values, Z_vectors = np.linalg.eig(Z)
print("The eigenvectors of X are\n\n", X_vectors,"\n")
print("And the eigenvalues of X are:", X_values)
# Python stores the eigenvectors as the columns of a matrix, so
# a corresponding eigenvalue-eigenvector pair is accessed
# like X_values[n], X_vectors[:,n] where n is 0 or 1.
# Compute the projectors for X, Y and Z:
proj_X_0 = X_vectors[:,0] * np.matrix.getH(X_vectors[:,0])
proj_X_1 = X_vectors[:,1] * np.matrix.getH(X_vectors[:,1])
proj_Y_0 = Y_vectors[:,0] * np.matrix.getH(Y_vectors[:,0])
proj_Y_1 = Y_vectors[:,1] * np.matrix.getH(Y_vectors[:,1])
proj_Z_0 = Z_vectors[:,0] * np.matrix.getH(Z_vectors[:,0])
proj_Z_1 = Z_vectors[:,1] * np.matrix.getH(Z_vectors[:,1])
print("The projector for +1 is:\n\n", proj_X_0,",\n")
print("and the projector for -1 is:\n\n", proj_X_1,".\n")
# The probabilities are computed as follows
p_X_0 = np.matrix.trace(dens * proj_X_0)
p_X_1 = np.matrix.trace(dens * proj_X_1)
p_Y_0 = np.matrix.trace(dens * proj_Y_0)
p_Y_1 = np.matrix.trace(dens * proj_Y_1)
p_Z_0 = np.matrix.trace(dens * proj_Z_0)
p_Z_1 = np.matrix.trace(dens * proj_Z_1)
print("The probability of getting +1 on X is: ", np.around(p_X_0,dp))
print("The probability of getting -1 on X is: ", np.around(p_X_1,dp))
print("The probability of getting +1 on Y is: ", np.around(p_Y_0,dp))
print("The probability of getting -1 on Y is: ", np.around(p_Y_1,dp))
print("The probability of getting +1 on Z is: ", np.around(p_Z_0,dp))
print("The probability of getting -1 on Z is: ", np.around(p_Z_1,dp))
