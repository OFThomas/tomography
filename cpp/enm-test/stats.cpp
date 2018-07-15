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
//
// Here, this is computed by finding
// the largest singular value, but
// I think it might be more complicated
// than that (as in, there are more
// variants of the norm to consider)
//
double distance_op(MatrixXc A, MatrixXc B) {
  Eigen::JacobiSVD<MatrixXc> svd(A-B, Eigen::ComputeThinU | Eigen::ComputeThinV);
  // The first singular value is the largest
  double distance = svd.singularValues()[0];
  return distance;
}

// Distance using the Frobenius norm
//
// This uses a library function from
// Eigen which computes the Frobenius
// (aka Hilbert-Schmidt) norm by
// default
//
double distance_trace(MatrixXc A, MatrixXc B){
  double distance = (A - B).norm();
  return distance;
}

// Fidelity distance
//
// The fidelity distance d is
// defined as follows:
//
//  F(A,B) = tr[ sqrt(sqrt(A) B sqrt(A)) ]
//  d(A.B) = arccos[F(A,B)]
//
// This is computed using a trick that
// involves diagonalising A to ease the
// sqrt function. This originated as a
// python optimisation -- it might not
// be necessary here.
double distance_fid(const MatrixXc A, const MatrixXc B) {

  // Obtain eigenvalues of A
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenA(A);
  if(eigenA.info() != Eigen::Success) abort();
  MatrixXc vectors_A = eigenA.eigenvectors();
  MatrixXc B_new = vectors_A.inverse() * B * vectors_A;
  // There's a bug hiding in here. The sqrt functions fails
  // if the input is a very small negative number. This needs to
  // be worked around somehow. I've fixed it by taking the
  // absolute value of the eigenvalue before sqrting. 
  MatrixXc A_new(2,2); A_new << std::sqrt(std::abs(eigenA.eigenvalues()[0])),0,
			 0, std::sqrt(std::abs(eigenA.eigenvalues()[1]));
  MatrixXc D = A_new * B_new * A_new;
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenD(D);
  if(eigenD.info() != Eigen::Success) abort();
  double fidelity = std::sqrt(std::abs(eigenD.eigenvalues()[0]))
    + std::sqrt(std::abs(eigenD.eigenvalues()[0])); 
  double distance = std::acos(fidelity);
  
  return distance;

}

// Fidelity distance 2
//
// This is a more direct implementation of the
// fidelity distance, directly from the definition.
//
double distance_fid_2(const MatrixXc A, const MatrixXc B) {
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenA(A);
  MatrixXc C = eigenA.operatorSqrt();
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenCBC(C * B * C);
  MatrixXc D = eigenCBC.operatorSqrt();
  double fidelity = std::real(D.trace());
  double distance = std::acos(fidelity);
  return distance;
}
