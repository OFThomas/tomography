/************************************************************
* Title: Testing the linear estimator
*
* Date created: 14th July 2018
*
* Language:    Cpp
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
*                 linear estimator
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

//#define DEBUG
//#define DEBUG_PRINT_RANDOM_DENSITY
//#define DEBUG_PRINT_PROBABILITY
//#define DEBUG_PRINT_DIAG
//#define DEBUG_PRINT_RANDOM
//#define DEBUG_PRINT_MEASUREMENTS
//#define DEBUG_PRINT_ESTIMATE
#define DEBUG_PRINT_ESTIMATES_EIGENVALUES

#include <iostream>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <random> 
#include <fstream>
#include <iomanip>
#include <chrono>
#include "simulation.h"
#include "estimation.h"
#include "stats.h"

int main() {

  auto start = std::chrono::steady_clock::now();
  
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

  //float av_distance[M][3];
  double non_physical[M];

  // Preliminaries: compute the projectors
  MatrixXc I(2,2); I << 1, 0, 0, 1;
  MatrixXc X(2,2); X << 0, 1, 1, 0;
  MatrixXc Y(2,2); Y << 0, std::complex<double>(0,-1),
		     std::complex<double>(0,1), 0;
  MatrixXc Z(2,2); Z << 1, 0, 0, -1;

  // Seed the random number generator
  srand(time(NULL));

  // Get an output file ready
  std::ofstream file;
  file.open("linear_test_1_cpp.dat");
  file << "Distances between estimated and "
       << "original density matrices using various distances."
       << "Go to the end of the file for the running time."
       << std::endl << std::endl;

  // Write the simulation parameters
  file << "Number of purity values tried = "
       << M << std::endl
       << "Number of density matrices per purity parameter = "
       << N << std::endl
       << "Total number of measurements for each of X, Y and Z = "
       << S << std::endl
       << std::endl;
  
  file << "PURITY, \tOPERATOR, \tTRACE, \t\tFIDELITY, \tNON PHYSICAL";
  // Set precision
  file << std::fixed << std::setprecision(5) << std::endl;
  
#ifdef DEBUG
#ifdef DEBUG_PRINT_IXYZ
  std::cout << "\n============== DEFINE X, Y AND Z ==============\n\n"; 
  std::cout << I << " This is I" << std::endl << std::endl;
  std::cout << X << " This is X" << std::endl << std::endl;
  std::cout << Y << " This is Y" << std::endl << std::endl;
  std::cout << Z << " This is Z" << std::endl << std::endl;
#endif
#endif

  MatrixXc vectors;
  MatrixXc proj_X[2];
  MatrixXc proj_Y[2];
  MatrixXc proj_Z[2];
  double outcomes_X[2];
  double outcomes_Y[2];
  double outcomes_Z[2];
  
  // Compute eigenvectors and eigenvalues of X, Y and Z
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenX(X);
  if(eigenX.info() != Eigen::Success) abort();
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenY(Y);
  if(eigenY.info() != Eigen::Success) abort();
  Eigen::SelfAdjointEigenSolver<MatrixXc> eigenZ(Z);
  if(eigenZ.info() != Eigen::Success) abort();

  // Compute all the projectors and store outcomes
  vectors = eigenX.eigenvectors();
  proj_X[0] = vectors.col(0) * vectors.col(0).adjoint();
  proj_X[1] = vectors.col(1) * vectors.col(1).adjoint();
  outcomes_X[0] = eigenX.eigenvalues()[0];
  outcomes_X[1] = eigenX.eigenvalues()[1];
  vectors = eigenY.eigenvectors();
  proj_Y[0] = vectors.col(0) * vectors.col(0).adjoint();
  proj_Y[1] = vectors.col(1) * vectors.col(1).adjoint();
  outcomes_Y[0] = eigenY.eigenvalues()[0];
  outcomes_Y[1] = eigenY.eigenvalues()[1];
  vectors = eigenZ.eigenvectors();
  proj_Z[0] = vectors.col(0) * vectors.col(0).adjoint();
  proj_Z[1] = vectors.col(1) * vectors.col(1).adjoint();
  outcomes_Z[0] = eigenZ.eigenvalues()[0];
  outcomes_Z[1] = eigenZ.eigenvalues()[1];
  
#ifdef DEBUG
#ifdef DEBUG_PRINT_EIGEN
  // Print all the eigenvalues and eigenvectors
  std::cout << "\n======== PRINT EIGENVALUES AND EIGENVECTORS ===========\n\n";
  std::cout << "The eigenvalues of X are\n"
	    << eigenX.eigenvalues()
	    << std::endl << std::endl;
  std::cout << "This matrix contains columns which are the eigenvectors of X\n"
	    << eigenX.eigenvectors()
	    << std::endl << std::endl;
  std::cout << "The eigenvalues of Y are\n"
	    << eigenY.eigenvalues()
	    << std::endl << std::endl;
  std::cout << "This matrix contains columns which are the eigenvectors of Z\n"
	    << eigenY.eigenvectors()
	    << std::endl << std::endl;
  std::cout << "The eigenvalues of Z are\n"
	    << eigenZ.eigenvalues(x)
	    << std::endl << std::endl;
  std::cout << "This matrix contains columns which are the eigenvectors of Z\n"
	    << eigenZ.eigenvectors()
	    << std::endl << std::endl;

  // Print all the projectors
  std::cout << "\n============== PRINT PROJECTORS ==============\n\n";
  std::cout << "The projector of X corresponding to outcome "
	    << eigenX.eigenvalues()[0] << " is\n"
	    << proj_X[0]
	    << std::endl << std::endl;
  std::cout << "The projector of X corresponding to outcome "
	    << eigenX.eigenvalues()[1] << " is\n"
	    << proj_X[1]
	    << std::endl << std::endl;
  std::cout << "The projector of Y corresponding to outcome "
	    << eigenY.eigenvalues()[0] << " is\n"
	    << proj_Y[0]
	    << std::endl << std::endl;
  std::cout << "The projector of Y corresponding to outcome "
	    << eigenY.eigenvalues()[1] << " is\n"
	    << proj_Y[1]
	    << std::endl << std::endl;
  std::cout << "The projector of Z corresponding to outcome "
	    << eigenZ.eigenvalues()[0] << " is\n"
	    << proj_Z[0]
	    << std::endl << std::endl;
  std::cout << "The projector of Z corresponding to outcome "
	    << eigenZ.eigenvalues()[1] << " is\n"
	    << proj_Z[1]
	    << std::endl << std::endl;
#endif
#endif
  
  // Define x -- put a loop here ------------------- LOOP for x between 0 and 1
  //
  // This loop runs through different values of the purity parameter x,
  // and tests the ability of the linear estimator in each case
  //

  // Variables to store the estimation error distances
  double dist_op[N];
  double dist_trace[N];
  double dist_fid[N];

  //int dp = 5; // Decimal places for printing
  double x = 0; // Purity parameter
  
  for(int k=0; k<M; k++) {
    // Temporary counter for non-physical estimates
    double non_physical_count = 0;				
    
    // Loop N times for each value of x ------ inner loop -- N trials for each x
    //
    // This loop generates N random density matrices for each fixed value of x
    // which used to simulate measurement data and run the estimator
    //
    for(int n=0; n<N; n++) {
      // Step 1: Prepare the density matrix
      //
      // The purity parameter x is picked between 0
      // and 1.
      //
      // Note: any time a numerical check is performed
      // and printed out, I've rounded the result to
      // make it more readable. Set the decimal places
      // to keep using the dp variable.
      //
      x = x_start + k * (x_end - x_start)/M;
      // Generate a density matrix with eigenvalues x and 1-x
      MatrixXc diag(2,2);
      diag << x, 0, 0, 1-x; 
      // Generate a random complex matrix
      MatrixXc cmat = MatrixXc::Random(2,2);
      // Obtain a random unitary matrix using Householder QR (see Eigen docs)
      Eigen::HouseholderQR<MatrixXc> qr(cmat);
      MatrixXc U = qr.householderQ(); // Get Q (the unitary matrix)
      MatrixXc dens = U * diag * U.adjoint();
      
#ifdef DEBUG
#ifdef DEBUG_PRINT_RANDOM
      std::cout << cmat << " This is a random complex matrix\n" << std::endl;
#endif
#ifdef DEBUG_PRINT_RANDOM_UNITARY
      std::cout << U << " This is a unitary matrix\n" << std::endl;
#endif
#ifdef DEBUG_PRINT_DIAG
      std::cout << diag << " This is the diagonal matrix\n" << std::endl;
#endif
#ifdef DEBUG_PRINT_RANDOM_DENSITY
      std::cout << dens << " This is the density matrix\n" << std::endl;
#endif
#endif
      // Step 2: Generate measurement data
      //
      // Generate data for X, Y and Z measurements
      //
      double X_data[S]; // To contain the (real) measurement values
      double Y_data[S];
      double Z_data[S];
      simulate(dens, proj_X, outcomes_X, S, X_data);
      simulate(dens, proj_Y, outcomes_Y, S, Y_data);
      simulate(dens, proj_Z, outcomes_Z, S, Z_data);

#ifdef DEBUG
#ifdef DEBUG_PRINT_MEASUREMENTS
      std::cout << "The simulated X values are:\n";
      for(int k=0; k<S; k++) std::cout << X_data[k] << ", ";
      std::cout << std::endl;
#endif
#endif
      
      // Step 3: Estimate density matrix
      //
      // Compute linear estimator
      //
      // Then tr(pI) is computed by requiring that
      // the density matrix be normalised
      //
      MatrixXc dens_est = linear_estimate_XYZ(X_data, Y_data, Z_data, S);

#ifdef DEBUG
#ifdef DEBUG_PRINT_ESTIMATE
      std::cout << dens_est << "This is the estimate"
		<< std::endl << std::endl;
#endif
#endif
      // Step 4: Compute and the distances
      //
      // Compute distances between the estimated
      // and true density matrix using the
      // different distance fuctions.
      //
      dist_op[n] = distance_op(dens, dens_est);
      dist_trace[n] = distance_trace(dens, dens_est);
      dist_fid[n] = distance_fid(dens, dens_est);
      // Count the number of non-physical matrices
      //
      Eigen::SelfAdjointEigenSolver<MatrixXc> eigenD(dens_est);
      if(eigenD.info() != Eigen::Success) abort();

#ifdef DEBUG
#ifdef DEBUG_PRINT_ESTIMATES_EIGENVALUES
      std::cout << "The estimated eigenvalues are "
		<< eigenD.eigenvalues()
		<< std::endl;
#endif
#endif
      if ((eigenD.eigenvalues()[0] < 0) || (eigenD.eigenvalues()[1] < 0)) {
	non_physical_count = non_physical_count + 1;
      }
      
    } // end of inner loop (fixed purity, random density matrices)

    // Step 5: Average the distances 
    //
    // Average the distances for each value of x. There are N density
    // matrices for each value of X, and consequently N distances to
    // compute. Therefore the mean is computed over N points.
    //
    double tmp_op(0), tmp_trace(0), tmp_fid(0);
    for(int k=0; k<N; k++) {
      tmp_op += dist_op[k];
      tmp_trace += dist_trace[k];
      tmp_fid += dist_fid[k];
    }
    double mean_op = tmp_op/N;
    double mean_trace = tmp_trace/N;
    double mean_fid = tmp_fid/N;
    //av_distance[k][0] = mean_op;
    //av_distance[k][1] = mean_trace;
    //av_distance[k][2] = mean_fid;
    non_physical[k] = non_physical_count/N;

    file << x << ",\t"
	 << mean_op << ",\t"
	 << mean_trace << ",\t"
	 << mean_fid << ",\t"
	 << non_physical[k]
	 << std::endl;
    
  } // end of outer loop (looping through different purities)

  // Get stop time
  auto end = std::chrono::steady_clock::now();
  auto time = end - start;

  auto dur = std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
  
  // Store the running time
  file << std::endl
       << "Total running time = "
       << dur << "ms" << std::endl;
    
  file.close();
}
