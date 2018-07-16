/************************************************************
* Title: Projectors and measurement operators
*
* Date created: 14th July 2018
*
* Language: Cpp
*
* Overview: This file contains various functions for
*           generating measurement projectors given
*           Hermitian operators
*
* Details: source
*
* Usage: Include the header file
*
************************************************************/

#include "proj.h"

// Generate projector
//
// Generator the projectors and outcomes for the Hermitian
// matrix A. The results are stored in two arrays. The
// order of the projectors in the proj_A array corresponds to
// the order of the outcomes in the outcomes_A array
//
// At the moment the function only works for 2x2 A.
//
int make_projector(MatrixXc A, MatrixXc proj_A[], double outcomes_A[]) {
  
  // Compute eigenvectors and eigenvalues of A
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenA(A);
  if(eigenA.info() != Eigen::Success) abort();
  
  // Store the projectors and outcomes
  MatrixXc vectors = eigenA.eigenvectors();
  proj_A[0] = vectors.col(0) * vectors.col(0).adjoint();
  proj_A[1] = vectors.col(1) * vectors.col(1).adjoint();
  outcomes_A[0] = eigenA.eigenvalues()[0];
  outcomes_A[1] = eigenA.eigenvalues()[1];

#ifdef DEBUG
#ifdef DEBUG_PRINT_PROJECTORS
  std::cout << std::endl
	    << "============================" << std::endl
	    <<"The projector of " << std::endl << std::endl
	    << A << std::endl << std::endl
	    << "corresponding to outcome "
	    << eigenA.eigenvalues()[0] << " is" << std::endl
	    << std::endl
	    << proj_A[0]
	    << std::endl << std::endl
	    << "and to outcome "
    	    << eigenA.eigenvalues()[1] << " is" << std::endl
	    << std::endl
	    << proj_A[1]
	    << std::endl;
#endif
#endif

  
  return 0;

}
