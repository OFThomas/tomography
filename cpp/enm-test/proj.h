/************************************************************
* Title: Projectors and measurement operators
*
* Date created: 15th July 2018
*
* Language: Cpp
*
* Overview: This file contains various functions for
*           generating measurement projectors given
*           Hermitian operators
*
* Details: Header
*
* Usage: Include the header file
*
************************************************************/

//#define DEBUG_PRINT_PROJECTORS

#include <iostream>
#include "Eigen/Dense"

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;

// Generate projector
//
// Generator the projectors and outcomes for the Hermitian
// matrix A. The results are stored in two arrays. The
// order of the projectors in the proj_A array corresponds to
// the order of the outcomes in the outcomes_A array
//
// At the moment the function only works for 2x2 A.
//
int make_projector(MatrixXc A, MatrixXc proj_A[], double outcomes_A[]);
