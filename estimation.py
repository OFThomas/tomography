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

# Function: linear_estimate_XYZ(meas_dat)
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
def linear_estimate_XYZ(meas_dat):
    I = np.matrix([[1,0],[0,1]])
    X = np.matrix([[0,1],[1,0]])
    Y = np.matrix([[0,-1j],[1j,0]])
    Z = np.matrix([[1,0],[0,-1]])
    mean_X = np.mean(meas_dat['data']['X'])
    mean_Y = np.mean(meas_dat['data']['Y'])
    mean_Z = np.mean(meas_dat['data']['Z'])
    dens_est = (mean_X * X + mean_Y * Y + mean_Z * Z + I)/2
    return dens_est
