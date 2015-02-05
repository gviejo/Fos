#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <math.h>
#include "pi_call.cpp"


using namespace std;

int main () {		
	int N = 52;		
	double fit [1] = {0.0};	
  	sferes_call(fit, N, "data_txt/B137_52/", 0.1, 0.1, 0.1);
  	std::cout << fit[0] << std::endl;

   	return 0;
}