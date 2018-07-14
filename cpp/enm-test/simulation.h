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
#define DEBUG_PRINT_PROBABILITY
//#define DEBUG_PRINT_DIAG
//#define DEBUG_PRINT_RANDOM
#define DEBUG_PRINT_MEASUREMENTS

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <random> 
#include "Eigen/Dense"

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;

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
	     double meas[],int S, double sim_dat[]);
