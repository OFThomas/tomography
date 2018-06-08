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

def estimate_rho(I, X, Y, Z, meas_dat, dp):
    mean_X = np.mean(meas_dat['X'])
    print("The mean of X is:", np.around(mean_X,dp))
    mean_Y = np.mean(meas_dat['Y'])
    print("The mean of Y is:", np.around(mean_Y,dp))
    mean_Z = np.mean(meas_dat['Z'])
    print("The mean of Z is:", np.around(mean_Z,dp))
    # The estimate for p is given by
    dens_est = (mean_X * X + mean_Y * Y + mean_Z * Z + I)/2
    return dens_est
