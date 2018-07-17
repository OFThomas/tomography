/************************************************************
* Title: ENM Test header file
*
* Date created: 15th July 2018
*
* Language: Cpp
*
* Overview: 
*
* Details: Header
*
* Usage: Include the header file
*
************************************************************/

//#define DEBUG_PRINT_RANDOM_DENSITY
//#define DEBUG_PRINT_PROBABILITY
//#define DEBUG_PRINT_DIAG
//#define DEBUG_PRINT_RANDOM
//#define DEBUG_PRINT_MEASUREMENTS
//#define DEBUG_PRINT_ESTIMATE
//#define DEBUG_PRINT_ESTIMATES_EIGENVALUES
//#define DEBUG_PRINT_EIGEN
#define DEBUG_PRINT_DISTANCE_AVERAGES

#define SHOW_PROGRESS

#include <iostream>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <random> 
#include <fstream>
#include <iomanip>
#include <chrono>
#include "Eigen/Dense"
#include "simulation.h"
#include "estimation.h"
#include "stats.h"
#include "proj.h"
#include "progress.h"

typedef Eigen::Matrix<std::complex<double>,
		      Eigen::Dynamic,
		      Eigen::Dynamic> MatrixXc;
