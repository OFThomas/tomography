/************************************************************
* Title: Summary statistics
*
* Date created: 19th June 2018
*
* Language:  Cpp
*
* Overview:    
*
* Details: source
*
* Usage:
*
************************************************************/

#include "stats.h"

// Distance using the operator norm
// Here, this is computed by finding
// the largest singular value, but
// I think it might be more complicated
// than that (as in, there are more
// variants of the norm to consider)
// 
float distance_op(MatrixXc A, MatrixXc B) {
  Eigen::JacobiSVD<MatrixXc> svd(A-B, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // The first singular value is the largest
  float distance = svd.singularValues()[0];
  return distance;
}
/*
// Distance using the Frobenius norm
//
float distance_trace(MatrixXc A, MatrixXc B){
  float distance = (A - B).norm();
  return distance;
}
    
def distance_fid(A,B):
    values_A,vectors_A = np.linalg.eig(A)
    B_new = np.matmul(np.matmul(np.linalg.inv(vectors_A), B), vectors_A)
    A_new = np.diag(np.sqrt(values_A))
    D = np.matmul(np.matmul(A_new, B_new), np.asmatrix(A_new))
    values_D = np.linalg.eigvals(D)
    #fidelity = np.sum(np.linalg.eigvals(D))
    fidelity = np.sum(np.sqrt(values_D))

    #C = sc.linalg.sqrtm(A) 
    #fidelity1 = np.matrix.trace(sc.linalg.sqrtm(C * B * C))
    #print(fidelity1-fidelity)
    #exit(1)
    distance = np.arccos(fidelity).real
    return distance

*/
