/************************************************************
* Title: File for printing progress
*
* Date created: 17th July 2018
*
* Language: Cpp
*
* Overview: 
*
* Details: header
*
* Usage: Include in the main program
*
************************************************************/

#include <iostream>
#include <string>
#include <sstream>
#include <chrono>

// Function: Show progress
//
// The function takes in a time since the
// start of the program and a proportion
// complete variable and prints a progress
// bar.
//
int show_progress(std::chrono::time_point<std::chrono::steady_clock> start,
		  double p);
