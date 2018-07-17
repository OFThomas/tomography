
#include "progress.h"

// Function: Show progress
//
// The function takes in a time since the
// start of the program and a proportion
// complete variable and prints a progress
// bar.
//
int show_progress(std::chrono::time_point<std::chrono::steady_clock> start,
		  double p) {
  
  std::string rem = "?";
  auto now = std::chrono::steady_clock::now();
  double t_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
  double t_s = t_ms/1000;
  std::ostringstream strs;
  strs << t_s/p - t_s;
  if(p != 0) rem = strs.str();
	      
  std::cout << "                      "
	    << "                      "
	    << "                    \r"
	    << "\u001b[0m\u001b[1m"
	    << "Percent done:\u001b[32m"
	    << (100 * p) 
	    << "%\u001b[0m\u001b[1m"
	    << "\t\tTotal time:\u001b[34m"
	    << t_s
	    << "s\u001b[0m\u001b[1m" 
	    << "\t\tEstimated time remaining:\u001b[35m"
	    << rem
	    << "s\u001b[0m\u001b[1m"
	    << "                 \r";
  if(p == 1) std::cout << std::endl;
  std::cout << std::flush;

  return 0;
}
