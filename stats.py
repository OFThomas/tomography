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
def distance_op(A,B):
    return distance = np.linalg.norm(A - B, 2)

def distance_trace(A,B):
    return distance = np.linalg.norm(A - B,'fro')

def distance_fid(A,B):
    fidelity = np.matrix.trace(sc.linalg.sqrtm(sc.linalg.sqrtm(A)
                    * B * sc.linalg.sqrtm(A)))        
    return distance = np.arccos(fidelity).real
