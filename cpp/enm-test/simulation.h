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
* Details: header
*
* Usage: Include in the main program
*
************************************************************/


#define DEBUG_PRINT_RANDOM_DENSITY
//#define DEBUG_PRINT_PROBABILITY
//#define DEBUG_PRINT_DIAG
//#define DEBUG_PRINT_RANDOM
#define DEBUG_PRINT_MEASUREMENTS

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <random> 
#include "Eigen/Dense"

#define _USE_MATH_DEFINES
#include <cmath>

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;

// Function: Generate random unitary
//
// The method is to parametrise the unitary
// group and then select the right distribution
// for the parameters. See '2009 Ozols - How to
// generate a random unitary matrix', page 5, for
// more details.
//
MatrixXc random_unitary(std::mt19937 & generator);

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
MatrixXc random_density(double x, std::mt19937 & generator);

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
int simulate(MatrixXc dens, const MatrixXc proj[],
	     const double meas[],int S, double sim_dat[],
	     std::mt19937 & generator);
