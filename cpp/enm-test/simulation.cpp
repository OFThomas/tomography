/************************************************************
* Title: Quantum state simulation
*
* Date created: 14th July 2018
*
* Language: Cpp
*
* Overview: This file contains various functions for
*           generating simulated states and measurement
*           data. 
*
* Details: source
*
* Usage: Include the header file
*
************************************************************/

#include "simulation.h"

// Function: random_density(x)
//
//  Generate a density matrix
//
//  The function takes a purity parameter x
//  which is between 0 and 1, which is one
//  of the eigenvalues of the density matrix
//  the other is 1-x.
//
//  It is better to generate a random
//  complex matrix and do a QR decomposition
//  than use the random unitary function.
//
MatrixXc random_density(double x) {

  // Generate a density matrix with eigenvalues x and 1-x
  MatrixXc diag(2,2);
  diag << x, 0, 0, 1-x; 
  // Generate a random complex matrix. The seed here is REALLY important!
  // In fact, this seeding business needs way more investigation.
  // It seems to completely change the resulting distributions, and
  // consequently the distribution of density matrix estimation errors
  // with purity value
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  srand(seed);      
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
  std::cout << std::endl
	    << "======================================"
	    << "======================================"
	    << std::endl
	    << "The randomly generated density matrix is"
	    << std::endl << std::endl
	    << dens
	    << std::endl << std::endl;
#endif
#endif

  return dens;
  
}
  
// Function: simulate()
//
//  Simulate measurement outcomes for
//  a set of projectors and associated
//  measurement outcomes.
//
//  The function generates S samples.
//
//  The array which will contain the
//  measurement simulations is passed
//  by reference.
//
int simulate(MatrixXc dens, const MatrixXc proj[],
	     const double meas[],int S, double sim_dat[]) {

  double p[2];
  for(int n=0; n<2; n++) {
    // Diagonalise it!!!
    p[n] = std::abs((dens * proj[n]).trace());
  }
  
#ifdef DEBUG
#ifdef DEBUG_PRINT_PROBABILITY
  std::cout << "Measurement probabilities:"
	    << std::endl << std::endl
	    << proj[0] << "   \tP = " << p[0]
	    << std::endl << std::endl
	    << proj[1] << "   \tP = " << p[1]
	    << std::endl << std::endl;
#endif
#endif

  // Just for fun, perform the following change:
  //
  // - unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  // - std::default_random_engine generator(seed);
  // + std::default_random_engine generator;
  //
  // and look at the change! Who knows what's going on with these random number
  // generators!
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::discrete_distribution<int> distribution {p[0], p[1]};
  // There must be a way of avoiding this if statement!
  for(int k=0; k<S; k++){
    int result = distribution(generator);
    if(result == 0) {
      sim_dat[k] = std::real(meas[0]);
    } else {
      sim_dat[k] = std::real(meas[1]);
    }
  }

#ifdef DEBUG
#ifdef DEBUG_PRINT_MEASUREMENTS
  std::cout << "Simulated measurements: ";
  for(int k=0; k<S; k++) std::cout << sim_dat[k] << ", ";
  std::cout << std::endl;
#endif
#endif
  
  
  return 0;
  
}
