#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <math.h>
#include "pi_call.cpp"

// extern "C" {
// 	int mvndst_(int *, double (*)[2], double (*)[2], int (*)[2], double *, int *,
//                 double *, double *, double *, double *, int *);
// }


using namespace std;

int main () {		
	int N = 52;		
	double fit [1] = {0.0};
	// int N = 2;
	// double lower [2] = {0.0,0.0};
	// double upper [2] = {1.0,1.0};
	// int infin [2] = {2,2};
	// double correl = 0.0;
	// int maxpts = 5000;
	// double abseps = 0.00005;
	// double releps = 0;
	// double error;
	// double value;
	// int inform;
	// mvndst_(&N, &lower, &upper, &infin, &correl, &maxpts, &abseps, &releps, &error, &value, &inform);
	// std::cout << value << " " << error << " " << inform << std::endl;
  	sferes_call(fit, N, "data_txt/B137_52/", 0.261945, 0.0227431, 0.1);
  	std::cout << fit[0] << std::endl;

   	return 0;
}