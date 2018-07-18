/************************************************************
* Title: Density matrix estimation
*
* Date created: 14th July 2018
*
* Language: Cpp
*
* Overview:    
*
* Details:
*
* Usage: to be imported
*
************************************************************/

#include "estimation.h"

// Function: linear_estimate_XYZ(X_data, Y_data, Z_data)
//
// The function takes in X, Y and Z
// measurement data and produces an
// estimate of the density matrix
// by computing the means of the X,
// Y and z data. The estimated density
// matrix is given by
//
// dens_est = 1/2 (mean_X // X +
//                 mean_Y // Y +
//                 mean_Z // Z +
//                 I) 
//
// This function returns a density
// matrix that satisfies the trace = 1
// constraint. This is achieved by
// leaving the multiple of I as a
// free parameter to be set after
// the other parameters have been
// estimated.
//
// The matrix is also Hermitian by
// construction. The only property
// not necessarily satisfied by the
// density matrix is positivity.
//
// S contains the number of samples
// in each of X, Y and Z
//
//
MatrixXc linear_estimate_XYZ(double X_data[],
			     double Y_data[],
			     double Z_data[],
			     int S) {
  MatrixXc I(2,2); I << 1, 0, 0, 1;
  MatrixXc X(2,2); X << 0, 1, 1, 0;
  MatrixXc Y(2,2); Y << 0, std::complex<double>(0,-1),
		     std::complex<double>(0,1), 0;
  MatrixXc Z(2,2); Z << 1, 0, 0, -1;

  // Compute means
  double mean_X = mean(X_data, S);
  double mean_Y = mean(Y_data, S);
  double mean_Z = mean(Z_data, S);
  
#ifdef DEBUG
#ifdef DEBUG_PRINT_MEANS
  std::cout << std::endl
	    << "X mean: " << mean_X << std::endl
	    << "Y mean: " << mean_Y << std::endl
	    << "Z mean: " << mean_Z << std::endl
	    << std::endl;
  
#endif
#endif
  
  // Reconstruct density matrix
  MatrixXc dens_est = (mean_X * X + mean_Y * Y + mean_Z * Z + I)/2;

#ifdef DEBUG
#ifdef DEBUG_PRINT_ESTIMATE
  std::cout << "The estimated density matrix is"
	    << std::endl << std::endl
	    << dens_est
	    << std::endl << std::endl;
#endif
#endif
  
  return dens_est;
}

// The objective function d is the distance between two matrices
// parametrised by X and Y. 
double d(const std::vector<double> & x, std::vector<double> & grad, void * f_data ) {
  double thing{*static_cast<double*>(f_data)};
  std::cout << thing;
  abort();
  return x[0] - thing;
}

// Function: enm_estimate_XYZ(X_data, Y_data, Z_data)
//
// This function estimates the density matrix using
// the extended norm minimisation method.
//
MatrixXc enm_estimate_XYZ(double X_data[],
			  double Y_data[],
			  double Z_data[]) {

  // Create an nlop object
  nlopt::opt opt(nlopt::LN_NELDERMEAD, 4); // 4 optimisation parameters
  // Set objective function
  int a{5};
  int * thing{&a};
  opt.set_min_objective(d, thing);
  std::vector<double> x{0,0,0,0};

  // Variable to contain the minimum of the function
  double minf;

  // Do the optimisation
  try {
    nlopt::result result = opt.optimize(x,minf);
    std::cout << "This is the output: " << minf << std::endl;
    std::cout << "The optimum x value is: "
	      << x[0] << ", " << x[1] << ", "
	      << x[2] << ", " << x[3] << std::endl;
    abort();
  }
  catch (std::exception & e) {
    std::cout << "nlopt failed: " << e.what() << std::endl;
    abort();
  }
    

    
  /*
  def d(x,rho_2):
            # x specifies rho = T T^
            rho_1 = np.matmul(np.array([[x[0], 0],[x[1]+1j*x[2], x[3]]]),np.array([[x[0], x[1]-1j*x[2]],[0, x[3]]]))
            distance = sts.distance_trace(rho_1, rho_2)
            return distance    

        rho_linear = linear_estimate_XYZ(X_data, Y_data, Z_data)
        result = sc.optimize.minimize(d, [1,1,1,1], args=rho_linear)
        x = result.x
        rho = np.matmul(np.array([[x[0], 0],[x[1]+1j*x[2], x[3]]]),np.array([[x[0], x[1]-1j*x[2]],[0, x[3]]]))
        return rho, rho_linear
  */

  
}
