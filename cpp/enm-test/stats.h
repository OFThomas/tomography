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
float distance_op(MatrixXc A, MatrixXc B);
