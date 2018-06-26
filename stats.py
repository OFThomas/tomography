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

def better_trace (A,N):
    trace = 0
    for n in range(0,N) : trace += A[n,n]
    return trace
    
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
    distance = np.linalg.norm(A - B, 2)
    return distance

def distance_trace(A,B):
    distance = np.linalg.norm(A - B,'fro')
    return distance
    
def distance_fid(A,B):
    values_A,vectors_A = np.linalg.eig(A)
    B_new = np.matmul(np.matmul(np.linalg.inv(vectors_A), B), vectors_A)
    A_new = np.diag(np.sqrt(values_A))
    D = np.asmatrix(A_new) * np.asmatrix(B_new) * np.asmatrix(A_new)
    values_D = np.linalg.eigvals(D)
    #fidelity = np.sum(np.linalg.eigvals(D))
    fidelity = np.sum(np.sqrt(values_D))

    C = sc.linalg.sqrtm(A) 
    fidelity1 = np.matrix.trace(sc.linalg.sqrtm(C * B * C))
    print(fidelity1-fidelity)
    exit(1)
    distance = np.arccos(fidelity).real
    return distance

# Use this one if you know the eigenvectors and eigenvalues of
# A already
def distance_fid_2(vectors_A,values_A,B):
    B_new = np.linalg.inv(np.matmul(np.matmul(vectors_A, B), vectors_A))
    A_new = np.diag(np.sqrt(values_A))
    D = np.asmatrix(A_new) * np.asmatrix(B_new) * np.asmatrix(A_new)
    values_D = np.linalg.eigvals(D)
    #fidelity = np.sum(np.linalg.eigvals(D))
    fidelity = np.sum(np.sqrt(values_D))

    #C = sc.linalg.sqrtm(A) 
    #fidelity1 = np.matrix.trace(sc.linalg.sqrtm(C * B * C))
    #print(fidelity1-fidelity)
    #exit(1)
    distance = np.arccos(fidelity).real
    return distance
