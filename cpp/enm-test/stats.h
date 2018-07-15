/************************************************************
* Title: Summary statistics
*
* Date created: 14th July 2018
*
* Language: Cpp
*
* Overview:    
*
* Details: Header file
*
* Usage: 
*
************************************************************/

#include "iostream"
#include "Eigen/Dense"
#include "Eigen/SVD"

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;

// Distance using the operator norm
// Here, this is computed by finding
// the largest singular value, but
// I think it might be more complicated
// than that (as in, there are more
// variants of the norm to consider)
// 
double distance_op(MatrixXc A, MatrixXc B);


// Distance using the Frobenius norm
//
// This uses a library function from
// Eigen which computes the Frobenius
// (aka Hilbert-Schmidt) norm by
// default
//
double distance_trace(MatrixXc A, MatrixXc B);

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
//
// distance_fid_2 implements the same
// operation in a different method
double distance_fid(const MatrixXc A, const MatrixXc B);
double distance_fid_2(const MatrixXc A, const MatrixXc B);

