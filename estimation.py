#############################################################
# Title: Density matrix estimation
#
# Date created: 7th June 2018
#
# Language: Python 3
#
# Overview:    
#
# Details:
#
# Usage: to be imported
#
#############################################################

import numpy as np
import scipy as sc

# Function: linear_estimate_XYZ(X_data, Y_data, Z_data)
#
# The function takes in X, Y and Z
# measurement data and produces an
# estimate of the density matrix
# by computing the means of the X,
# Y and z data. The estimated density
# matrix is given by
#
# dens_est = 1/2 (mean_X * X +
#                 mean_Y * Y +
#                 mean_Z * Z +
#                 I) 
#
# This function returns a density
# matrix that satisfies the trace = 1
# constraint. This is achieved by
# leaving the multiple of I as a
# free parameter to be set after
# the other parameters have been
# estimated.
#
# The matrix is also Hermitian by
# construction. The only property
# not necessarily satisfied by the
# density matrix is positivity.
#
def linear_estimate_XYZ(X_data, Y_data, Z_data):
    I = np.matrix([[1,0],[0,1]])
    X = np.matrix([[0,1],[1,0]])
    Y = np.matrix([[0,-1j],[1j,0]])
    Z = np.matrix([[1,0],[0,-1]])
    mean_X = np.mean(X_data)
    mean_Y = np.mean(Y_data)
    mean_Z = np.mean(Z_data)
    dens_est = (mean_X * X + mean_Y * Y + mean_Z * Z + I)/2
    return dens_est

# Function: maximum_likelihood_XYZ(X_data, Y_data, Z_data)
#
# This function estimates the density matrix using
# X, Y and Z data using the maximum likelihood method
#
# The likelihood of the density matrix rho given
# N measurements outcomes P_n is
#
#              _N_
#             |   |
#    L[rho] = |   | Tr[rho * P_n]  
#             n = 0
#
# There are six projectors for X, Y and Z;
# one for each outcome +1 or -1. If X, Y and
# Z are each measured S times, then N = 3S and
# there are N terms in the product. The first
# S terms correspond to X, the second to Y and
# the last to Z. 
#
# The proj array contains the X, Y and Z projectors
# in the order X, Y and Z, with + before -. 
#
# The parameter space is specified as follows:
# Let rho = T^ T where T is a lower triangular
# complex matrix with real elements on the
# diagonal. This guarantees that rho is
# Hermitian and positive. Four parameters are
# necessary to specify T in the one qubit case.
#
# The trace = 1 condition is fulfilled using
# the method of lagrange multipliers.
#
def maximum_likelihood_XYZ(X_data, Y_data, Z_data, proj):
    def L(a,b,c,d,lag):
        # a,b,c,d specify rho = T^ T
        rho = [[a, b-1j*c],[0, d]] * [[a, 0],[b+1j*c, d]]
        return tmp
        
    result = sc.optimize.minimize(L, 1)
    return result


