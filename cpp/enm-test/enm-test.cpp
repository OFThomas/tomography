/************************************************************
* Title: Testing extended norm minimisation
*
* Date created: 11th June 2018
*
* Language:    C++
*
* Overview:    The program tests the extended norm
*              minimisation algorithm in the case
*              of an informationally complete set of
*              measurements 
*
* Details:     The script performs the following steps:
*
*              1) Set x in [0,1] as the purity parameter;
*                 pick a random unitary U; and compute
*
*
*                          -       -
*                         | x     0 |
*                 p = U * |         | * U^
*                         | 0   1-x |
*                          -       -
*
*                 The closer x is to 1 or 0, the more
*                 pure is the state.
*
*              2) Fix a set of measurement operators,
*                 and generate sample data from the
*                 distribution of outcomes implied by
*                 the state p.
*
*              3) Compute the density matrix using the
*                 extended norm minimsation.
*
*              4) Compute the distance between the
*                 estimate and the true density matrix
*                 for each of the different distances
*                 (operator norm, Hilbert-Schmidt norm,
*                 fidelity).
*
*              5) Perform averages of each distance at each
*                 value of x. Plot the average distances
*                 as a function of x for all the values of
*                 x.
*
*              6) Repeat steps 1-4 for different values
*                 of x and repeat each value of x a large
*                 number of times. Store all the distances.
*
* Compilation: make 
*
*************************************************************/

#define DEBUG

#include <iostream>
#include <complex>
#include "Eigen/Dense"

int main() {

  using Eigen::MatrixXd;
  typedef Eigen::Matrix<std::complex<double>,
  			Eigen::Dynamic,
  			Eigen::Dynamic> MatrixXc;
  
  // ======= Test parameter ===============================
  int M = 2000;  // Number of purity parameters x to try
  double x_start = 0; // Specify purity parameter x range
  double x_end = 1;
  int N = 500;  // Number of random density matrices per x value
  int S = 500;  // Number of samples of each measurement to
  // simulate for each density matrix 
  // ======================================================

  float av_distance[M][3];
  float non_physical[M];

  // Preliminaries: compute the projectors
  MatrixXc I(2,2); I << 1, 0, 0, 1;
  MatrixXc X(2,2); X << 0, 1, 1, 0;
  MatrixXc Y(2,2); Y << 0, (0,-1), (0,1), 0;
  MatrixXc Z(2,2); Z << 1, 0, 0, -1;

  
#ifdef DEBUG
  std::cout << I << " This is I" << std::endl << std::endl;
  std::cout << X << " This is X" << std::endl << std::endl;
  std::cout << Y << " This is Y" << std::endl << std::endl;
  std::cout << Z << " This is Z" << std::endl << std::endl;
#endif

  // values_X, vectors_X = np.linalg.eig(X)
  // proj_X = np.zeros([2,2,2])
  // proj_X[0,:,:] = np.matmul(vectors_X[:,0], np.matrix.getH(vectors_X[:,0]))
  // proj_X[1,:,:] = np.matmul(vectors_X[:,1], np.matrix.getH(vectors_X[:,1]))
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigensolver(X);
  if(eigensolver.info() != Eigen::Success) abort();
#ifdef DEBUG
  std::cout << "The eigenvalues of X are\n"
	    << eigensolver.eigenvalues()
	    << std::endl;
  std::cout << "This matrix contains columns which are the eigenvectors of X\n"
	    << eigensolver.eigenvectors()
	    << std::endl;
#endif
  
}
