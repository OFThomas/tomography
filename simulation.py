#############################################################
# Title: Quantum state simulation
#
# Date created: 7th June 2018
#
# Language: Python 3
#
# Overview: This file contains various functions for
#           generating simulated states and measurement
#           data. 
#
# Details:
#
# Usage: to be imported
#
#############################################################

import numpy as np
from scipy.stats import unitary_group as ug
from input import *

# Function: random_density(x)
#
#   Generate a density matrix
#
#   The function takes a purity parameter x
#   which is between 0 and 1, which is one
#   of the eigenvalues of the density matrix
#   the other is 1-x.
#
def random_density(x):
    # Check that the purity parameter is valid
    if x > 1 or x < 0:
        print("Bug: invalid input argument x:",x)
        print("Exiting")
        exit(1)
    realMat = np.random.random((2,2))
    U = np.asmatrix(ug.rvs(2)) # Random unitary
    # Check that U is actually unitary
    U_dag = np.matrix.getH(U)
    if(U_dag * U - np.asmatrix([[1,0],[0,1]]) > 1e-5).any(): 
        print("Bug: failed to generate unitary matrix")
        print("Exiting")
        exit(1)
    # Compute the density matrix
    diag = np.matrix([[x,0],[0,1-x]])
    dens = U * diag * U_dag
    # Check the density matrix
    if abs(np.trace(dens) - 1) > 1e-10:
        print("Bug: failed to generate a density matrix with trace 1")
        print("Exiting")
        exit(1)
    if np.linalg.eig(dens)[0][0] < -1e-5 or np.linalg.eig(dens)[0][1] < -1e-5:
        print("Bug: failed to generate a positve density matrix")
        print("Exiting")
        exit(1)
    return dens



# Function: density(dp)
#
#   Generate a density matrix. The functions
#   provides the user with the option to 
#   enter the matrix manually, or it generates
#   the matrix at random. 
#
#   The function takes a purity parameter x
#   which is between 0 and 1.
#
def density(x,dp):
    print("===== Density matrix generation =====\n")
    print("A density matrix is required to test")
    print("the state tomography proceedure. It")
    print("can be input manually or generated")
    print("randomly.")
    response = yes_or_no("Do you want to enter a density matrix manually? ")
    if response == 0:
        # Request matrix dimension
        while 1==1:
            dim = get_user_value("Specify the matrix dimension", "integer")
            if dim >= 2: 
                if bin(dim).count("1") == 1: break
                else: print("Matrix dimension should be power of 2")  
            else: print("Matrix dimension must be at least 2")
        # Populate matrix
        dens = np.asmatrix(np.zeros((dim,dim), dtype=complex))
        while 1==1:
            for m in np.ndindex((dim,dim)):
                dens[m] = get_user_value("Input element " + str(m), "complex")
            # Print and check the density matrix
            response = yes_or_no("Is this the correct density matrix?\n\n"+str(dens)+"\n\n")
            if response == 0: return dens
            else: print("Enter the density matrix again.")
    else: 
        # Generate a random matrix
        print("\nThe random density matrix will be generated")
        print("by generating random eignevalues x and 1-x")
        print("and then conjugating the resulting diagonal")
        print("matrix by a random unitary matrix U\n")
        dens = random_density(x)
        return dens

# Function: density(dp)
#
#   Generate a density matrix. The functions
#   provides the user with the option to 
#   enter the matrix manually, or it generates
#   the matrix at random. 
#
#   The function takes a purity parameter x
#   which is between 0 and 1.
#
#   This assumes that the measurement operators
#   have distinct eigenvalues (no algebraic
#   multiplicities).
#
#   The function returns a dictionary with
#   two main branches: 'data' and
#   'probabilities'. The first contains a
#   dictionary of simulated data indexed
#   by measurement. The second contains a
#   dictionary of probabilities indexed
#   first by measurement and then by
#   measurement outcome.
#
def simulate(dens,meas_ops,samples):
    values, vectors, eigen, proj, sim_dat = {},{},{},{},{}
    sim_dat['probabilities'], sim_dat['data'] = {},{}
    # Compute the eigenvectors and eigenvalues of X, Y and Z
    for measurement in meas_ops:
        proj[measurement],eigen[measurement] = {},{}
        sim_dat['probabilities'][measurement] = {}
        values[measurement], vectors[measurement] = np.linalg.eig(np.asmatrix(meas_ops[measurement]))
        for k in range(0,values[measurement].size):
            eigen[measurement][values[measurement][k]] = vectors[measurement][:,k]
        # Compute the projectors for X, Y and Z:
        p,v = [],[]
        for outcome in eigen[measurement]:
            sim_dat['data'][measurement] = {}
            proj[measurement][outcome] = eigen[measurement][outcome] * np.matrix.getH(eigen[measurement][outcome])
            sim_dat['probabilities'][measurement][outcome] = np.trace(dens * proj[measurement][outcome]).real 
            # Generate the measurement data
            p.append(sim_dat['probabilities'][measurement][outcome])
            v.append(complex(outcome))
        sim_dat['data'][measurement] = np.random.choice(v,samples,p=p)
    return sim_dat
