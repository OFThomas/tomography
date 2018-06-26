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
#   It is better to generate a random
#   complex matrix and do a QR decomposition
#   than use the random unitary function.
#
def random_density(x):
    cmat = np.random.random([2,2]) + 1j * np.random.random([2,2])
    U = np.linalg.qr(cmat)[0]
    #U = ug.rvs(2) # Random unitary
    U_dag = np.matrix.getH(np.asmatrix(U))
    diag = [[x,0],[0,1-x]]
    dens = np.matmul(np.matmul(U, diag), U_dag)
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

def better_trace (A,N):
    trace = 0
    for n in range(0,N) : trace += A[n,n]
    return trace

def better_trace_2 (A):
    trace = A[0,0] + A[1,1]
    return trace

    
# Function: simulate()
#
#   Simulate measurement outcomes for
#   a set of projectors and associated
#   measurement outcomes. The proj
#   variable contains the projectors and
#   the meas variable contains the
#   outcomes. Both variables are arrays
#   ordered so that the ith outcome in
#   meas corresponds to the ith projector
#   in proj.
#
#   The first index of proj indexes the
#   projectors, and the second two
#   index the matrix.
#
#   The function generates S samples.
#
def simulate(dens,proj,meas,S):
    N = meas.size
    sim_dat = np.zeros(S)
    p = np.zeros(N)
    for n in range(N):
        ## Diagonalise it!!!
        p[n] = better_trace_2(np.matmul(dens, proj[n,:,:])).real
    #sim_dat = np.random.choice(meas,S,p=p)
    tmp = np.random.rand(S)
    # Need to update this for more measurement outcomes!
    #for k in range(S):
        #if tmp[k] < p[0] : sim_dat[k] = meas[0]
        #else : sim_dat[k] = meas[1]
    sim_dat[tmp <= p[0]] = meas[0]
    sim_dat[tmp > p[0]] = meas[1]    
    return sim_dat
