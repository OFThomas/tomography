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

// Function: Generate random unitary
//
// The method is to parametrise the unitary
// group and then select the right distribution
// for the parameters. See '2009 Ozols - How to
// generate a random unitary matrix', page 5, for
// more details.
//
MatrixXc random_unitary(std::mt19937 & generator) {

  // Pick alpha, phi and chi uniformly in the
  // interval [0, 2*pi]
  std::uniform_real_distribution<double> uniform_1(0.0, 2 * M_PI);
  double alpha = uniform_1(generator);
  double psi = uniform_1(generator);
  double chi = uniform_1(generator);
  // Pick xi uniformly from [0,1]
  std::uniform_real_distribution<double> uniform_2(0.0, 1.0);
  double xi = uniform_2(generator);

  // Compue derived quantities
  double phi = std::asin(std::sqrt(xi));
  std::complex<double> glob =  std::exp(std::complex<double>(0,alpha));
  
  // Compute matrix elements
  std::complex<double> a = std::exp(std::complex<double>(0,psi)) * std::cos(phi); 
  std::complex<double> b = std::exp(std::complex<double>(0,chi)) * std::sin(phi);

  // Write the matrix
  MatrixXc U(2,2); U << a, b, -std::conj(b), std::conj(a);
  // Global phase
  U = glob * U;

#ifdef DEBUG
#ifdef DEBUG_PRINT_UNITARY
  std::cout << U << "This is it" << std::endl
	    << U * U.adjoint() << "What is this?"
	    << std::endl;
#endif
#endif
  return U;
}

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
MatrixXc random_density(double x, std::mt19937 & generator) {

  // Generate a density matrix with eigenvalues x and 1-x
  MatrixXc diag(2,2);
  diag << x, 0, 0, 1-x; 
  // Generate a random complex matrix. The seed here is REALLY important!
  // In fact, this seeding business needs way more investigation.
  // It seems to completely change the resulting distributions, and
  // consequently the distribution of density matrix estimation errors
  // with purity value. The seed is generated at the start of the
  // program
  //
  // There's a subtelty here, which is that it is important to use
  // the Gram-Schmidt method for performing the QR decomposition.
  // (See '2009 Ozols - How to generate a random unitary matrix', page 6.)
  //
  //std::complex<double> a = std::complex<double>(1,1); 
  //MatrixXc correction(2,2); correction << a,a,a,a;
  //MatrixXc cmat = MatrixXc::Random(2,2);
  // Obtain a random unitary matrix using Householder QR (see Eigen docs)
  //Eigen::HouseholderQR<MatrixXc> qr(cmat);
  //MatrixXc U = qr.householderQ(); // Get Q (the unitary matrix)
  MatrixXc U = random_unitary(generator);
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
//  Note that the random number
//  generator in the last argument
//  must be passed by reference.
//  Otherwise a copy gets made and
//  the random numbers from different
//  calls will end up dependent on each
//  other
//
// It appears that there is no difference
// between using default_random_engine and
// mt19937 when it comes to the general
// look of the output distributions.
//
int simulate(MatrixXc dens, const MatrixXc proj[],
	     const double meas[],int S, double sim_dat[],
	     std::mt19937 & generator) {

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
  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
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
