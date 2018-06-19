#############################################################
# Title: Summary statistics
#
# Date created: 19th June 2018
#
# Language:    Python 3
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

# Distance between density matrices
#
# Compute the distance between A and
# B using some distance specified by
# type.
#
# A metric d is always obtained from
# a norm using
#
#     d(A,B) = ||A - B||
#
# 1) If type is 'operator', compute the
#    distance based on the operator norm:
#
#    The operator norm of A is
#
#        ||A|| = supp |Ax|
#               |x|=1
#
# 2) If type is 'trace', compute the
#    distance based on the Hilbert-Schmidt
#    norm:
#
#    The Hilbert-Schmidt norm of A is:
#
#        ||A|| = sqrt[ tr(A^ A) ]
#
# 3) If type is 'fidelity', compute the
#    fidelity distance is defined as
#    follows
#
#        F(A,B) = tr[ sqrt(sqrt(A) B sqrt(A)) ]
#        d(A.B) = arccos[F(A,B)]
#
def distance(A,B,type):
    if type == 'operator':
        distance = np.linalg.norm(A - B, 2)
    elif type == 'trace':
        distance = np.linalg.norm(A - B,'fro')
    elif type == 'fidelity':
        fidelity = np.matrix.trace(sc.linalg.sqrtm(sc.linalg.sqrtm(A)
                    * B * sc.linalg.sqrtm(A)))        
        distance = np.arccos(fidelity).real
    return distance

'''
variances = {}
# Get the purity of the density matrix estimate
eigenvalues, eigenvectors = np.linalg.eig(dens_est)
# Assume eigenvalues are real
print("The eigenvalues of the estimates density matrix are:\n\n\t",
      np.around(eigenvalues.real[0],dp), "and",
      np.around(eigenvalues.real[1],dp))
print("\nThe purity parameter of the estimated density matrix is",
      np.maximum(np.around(eigenvalues.real[0],dp),
      np.around(eigenvalues.real[1],dp)))
print("\nThe variances in the simulated data are:\n")
for key in meas_dat :
    variances[key] = np.var(meas_dat[key])
    print("\tThe variance in the",key,"samples was",variances[key])
'''
