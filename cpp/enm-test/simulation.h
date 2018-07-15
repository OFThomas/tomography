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

//#define DEBUG
//#define DEBUG_PRINT_RANDOM_DENSITY
//#define DEBUG_PRINT_PROBABILITY
//#define DEBUG_PRINT_DIAG
//#define DEBUG_PRINT_RANDOM
//#define DEBUG_PRINT_MEASUREMENTS

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <random> 
#include "Eigen/Dense"

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;

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
MatrixXc random_density(double x);

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
	     const double meas[],int S, double sim_dat[]);
