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
import stats as sts

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


class ML:

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
    # the method of lagrange multipliers. If
    #
    #          _       _
    #         | a     0 |
    #     T = |         |
    #         | b+jc  d |
    #          ^       ^
    #
    @staticmethod
    def XYZ(X_data, Y_data, Z_data, proj):
    
        def L(x,args):
            # x specifies rho = T T^
            rho = np.matmul(np.array([[x[0], 0],[x[1]+1j*x[2], x[3]]]),np.array([[x[0], x[1]-1j*x[2]],[0, x[3]]]))
            tmp = 1
            M = args[:,0,0].size
            for n in range(0,M):
                tmp = tmp * np.real(np.trace(np.matmul(rho, args[n,:,:])))
            tmp = tmp - x[4] * np.real(np.trace(rho)) # Lagrange multiplier
            return -tmp

        def L_1(x,args):
            # x specifies rho = T T^
            rho = np.matmul(np.array([[x[0], 0],[x[1]+1j*x[2], x[3]]]),np.array([[x[0], x[1]-1j*x[2]],[0, x[3]]]))
            print(rho)
            tmp = 0
            M = args[:,0,0].size
            print(args)
            exit(1)
            for n in range(0,M):
                tmp = tmp + np.log(np.real(np.trace(np.matmul(rho, args[n,:,:]))))
                #print(np.trace(np.matmul(rho, args[n,:,:])))
                #exit(1)
            tmp = tmp - x[4] * np.real(np.trace(rho)) # Lagrange multiplier
            print("tmp:",-tmp)
            return -tmp

        def L_2(x,args):
            # x specifies rho = T T^
            rho = np.matmul(np.array([[x[0], 0],[x[1]+1j*x[2], x[3]]]),np.array([[x[0], x[1]-1j*x[2]],[0, x[3]]]))
            M = args[:,0,0].size
            trace = np.zeros(M)
            for n in range(0,M):
                trace[n] = np.log(np.real(np.trace(np.matmul(rho, args[n,:,:]))))
            total = -trace.sum() + x[4] * np.real(np.trace(rho)) # Lagrange multiplier
            print("tot:",-total)
            return total


        data = np.concatenate((X_data, Y_data, Z_data), axis=0)
        #print("X",X_data)
        #print("Y",Y_data)
        #print("Z",Z_data)
        S = X_data.size
        N = 3*S
        args = np.zeros([N,2,2],dtype='complex')
        for n in range(S):
            if X_data[n] == +1 : args[n,:,:] = proj[0]
            else : args[n,:,:] = proj[1]
            if Y_data[n] == +1 : args[n+S,:,:] = proj[2]
            else : args[n+S,:,:] = proj[3]
            if Z_data[n] == +1 : args[n+2*S,:,:] = proj[4]
            else : args[n+2*S,:,:] = proj[5]

        #print(args)
        # fun1 = lambda var : var[0]**2 + var[1]**2 + var[2]**2 + var[3]**2 - 1
        # cons = ({'type':'eq', 'fun':fun1},
        #         {'type':'eq', 'fun':lambda x : (x[0]**2)*(x[3]**2) - 0.25},
        #         {'type':'ineq', 'fun':lambda x : x[0]},
        #         {'type':'ineq', 'fun':lambda x : x[3]})

        # cons1 = ({'type':'ineq', 'fun':lambda x : x[0]},
        #          {'type':'ineq', 'fun':lambda x : x[3]})

        #result = sc.optimize.minimize(L_2, [1,0,0,0,0], args=args, constraints=cons1, method='TNC')
        #result = sc.optimize.minimize(L_2, [1,0,0,0,0], args=args)
        x = sc.optimize.fmin(L_2, [1,0,0,0,0], args=(args,))
        exit(1)
        #print(result)
        #exit(1)
        #x = result.xopt
        rho = np.matmul(np.array([[x[0], 0],[x[1]+1j*x[2], x[3]]]),np.array([[x[0], x[1]-1j*x[2]],[0, x[3]]]))
        return rho



# Function: enm_XYZ(X_data, Y_data, Z_data)
#
# This function estimates the density matrix using
# the extended norm minimisation method.
#
# 
#
class ENM():

    @staticmethod
    def enm_XYZ(X_data, Y_data, Z_data):

        def d(x,rho_2):
            # x specifies rho = T T^
            rho_1 = np.matmul(np.array([[x[0], 0],[x[1]+1j*x[2], x[3]]]),np.array([[x[0], x[1]-1j*x[2]],[0, x[3]]]))
            distance = sts.distance_trace(rho_1, rho_2)
            return distance    

        rho_linear = linear_estimate_XYZ(X_data, Y_data, Z_data)
        result = sc.optimize.minimize(d, [1,1,1,1], args=rho_linear)
        x = result.x
        rho = np.matmul(np.array([[x[0], 0],[x[1]+1j*x[2], x[3]]]),np.array([[x[0], x[1]-1j*x[2]],[0, x[3]]]))
        return rho, rho_linear
