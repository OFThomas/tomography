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
int simulate(MatrixXc dens, MatrixXc proj[],
	     double meas[],int S, double sim_dat[]) {

  double p[2];
  for(int n=0; n<2; n++) {
    // Diagonalise it!!!
    p[n] = std::abs((dens * proj[n]).trace());
  }
  
#ifdef DEBUG
#ifdef DEBUG_PRINT_PROBABILITY
  std::cout << "Probability distribution for X:\n"
	    << p[0] << ", "
	    << p[1] << std::endl << std::endl;
#endif
#endif
  
  std::default_random_engine generator;
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
  
  return 0;
  
}
