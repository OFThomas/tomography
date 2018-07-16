/************************************************************
* Title: Density matrix estimation
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

#define DEBUG_PRINT_MEANS
#define DEBUG_PRINT_ESTIMATE

#include <iostream>
#include "stats.h"
#include "Eigen/Dense"

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;

// Function: linear_estimate_XYZ(X_data, Y_data, Z_data)
//
// The function takes in X, Y and Z
// measurement data and produces an
// estimate of the density matrix
// by computing the means of the X,
// Y and z data. The estimated density
// matrix is given by
//
// dens_est = 1/2 (mean_X // X +
//                 mean_Y // Y +
//                 mean_Z // Z +
//                 I) 
//
// This function returns a density
// matrix that satisfies the trace = 1
// constraint. This is achieved by
// leaving the multiple of I as a
// free parameter to be set after
// the other parameters have been
// estimated.
//
// The matrix is also Hermitian by
// construction. The only property
// not necessarily satisfied by the
// density matrix is positivity.
///
MatrixXc linear_estimate_XYZ(double X_data[],
			     double Y_data[],
			     double Z_data[],
			     int S);
