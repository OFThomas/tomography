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

#include "enm-test.h"

int main() {

  // Start the clock!
  std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();

  
  // Step 0: Seed the random number generators
  //
  // The random engine is initialised outside the main program
  // This is so that the engine is only initialised once,
  // as opposed to every time the simulate function is
  // called. That way each set of measurements is independent
  // of the others. In fact, it is critical that it is passed
  // by reference -- otherwise the object gets copied and
  // you get the same problem. If the generator is passed
  // by reference, the same object is used by all the
  // random number generators.
  //
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  srand(seed);
  std::mt19937 gen;
  gen.seed(seed);

  // ======= Test parameter ===============================
  int M = 2000;  // Number of purity parameters x to try
  double x_start = 0; // Specify purity parameter x range
  double x_end = 1;
  int N = 500;  // Number of random density matrices per x value
  int S = 500;  // Number of samples of each measurement to
  // simulate for each density matrix 
  // ======================================================

  double non_physical[M];

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

  // Preliminaries: define measurement operators
  MatrixXc I(2,2); I << 1, 0, 0, 1;
  MatrixXc X(2,2); X << 0, 1, 1, 0;
  MatrixXc Y(2,2); Y << 0, std::complex<double>(0,-1),
		     std::complex<double>(0,1), 0;
  MatrixXc Z(2,2); Z << 1, 0, 0, -1;

  // Compute the projectors
  MatrixXc proj_X[2];
  MatrixXc proj_Y[2];
  MatrixXc proj_Z[2];
  double outcomes_X[2];
  double outcomes_Y[2];
  double outcomes_Z[2];
  
  make_projector(X, proj_X, outcomes_X);
  make_projector(Y, proj_Y, outcomes_Y);
  make_projector(Z, proj_Z, outcomes_Z);
  
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
#ifdef DEBUG
    std::cout << "======================= "
	      << "START OF AN X LOOP"
	      << "======================="
	      << std::endl;
#endif
    // Temporary counter for non-physical estimates
    double non_physical_count = 0;
    
    // Loop N times for each value of x ------ inner loop -- N trials for each x
    //
    // This loop generates N random density matrices for each fixed value of x
    // which used to simulate measurement data and run the estimator
    //
    for(int n=0; n<N; n++) {
      #ifdef DEBUG
      std::cout << std::endl
		<< "++++++++ "
		<< "FIXED DENSITY MATRIX "
		<< "++++++++"
		<< std::endl << std::endl;
#endif
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
      MatrixXc dens = random_density(x, gen);
      //MatrixXc dens(2,2); dens << 1,0,0,0;
      
      // Step 2: Generate measurement data
      //
      // Generate data for X, Y and Z measurements. 
      //
      double X_data[S]; // To contain the (real) measurement values
      double Y_data[S];
      double Z_data[S];
      simulate(dens, proj_X, outcomes_X, S, X_data, gen);
      simulate(dens, proj_Y, outcomes_Y, S, Y_data, gen);
      simulate(dens, proj_Z, outcomes_Z, S, Z_data, gen);
      
      // Step 3: Estimate density matrix
      //
      // Compute linear estimator
      //
      // Then tr(pI) is computed by requiring that
      // the density matrix be normalised
      //
      //MatrixXc dens_est = linear_estimate_XYZ(X_data, Y_data, Z_data, S);
      MatrixXc dens_est = enm_estimate_XYZ(X_data, Y_data, Z_data, S);

      // Step 4: Compute and the distances
      //
      // Compute distances between the estimated
      // and true density matrix using the
      // different distance fuctions.
      //
      dist_op[n] = distance_op(dens_est, dens);
      dist_trace[n] = distance_trace(dens_est, dens);
      dist_fid[n] = distance_fid_2(dens_est, dens);

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
    double mean_op = mean(dist_op, N);
    double mean_trace = mean(dist_trace, N);
    double mean_fid = mean(dist_fid, N);
    non_physical[k] = non_physical_count/N;

#ifdef SHOW_PROGRESS
    // Show progress
    double p = static_cast<double>(k+1)/M;
    show_progress(start,p);
#endif
    
    
#ifdef DEBUG
#ifdef DEBUG_PRINT_DISTANCE_AVERAGES
    std::cout << std::endl
	      << "The average operator distance is: " << mean_op << std::endl    
	      << "The average trace distance is: " << mean_trace << std::endl
	      << "The average fidelity distance is: " << mean_fid << std::endl;
#endif
#endif
    
    // Step 6: Write the results to a file
    //
    file << x << ",\t"
	 << mean_op << ",\t"
	 << mean_trace << ",\t"
	 << mean_fid << ",\t"
	 << non_physical[k]
	 << std::endl;

#ifdef DEBUG
    std::cout << std::endl
	      << "++++++++ "
	      << "END OF INNER LOOP "
	      << "++++++++"
	      << std::endl << std::endl;
#endif
    
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
