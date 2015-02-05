#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <math.h>
#include "pi_call.cpp"
// #include "normal.cpp"

using namespace std;

int main () {		
	int N = 52;		
	double fit [1] = {0.0};	
  	sferes_call(fit, N, "data_txt/B137_52/", 0.5, 0.1, 0.1);  	
  	std::cout << fit[0] << std::endl;


	// int n_case = 30;
	// double grid [900] [2];
	// double ystart = -3.0;
	// double p_goal [900];
	// for (int i=0;i<n_case;i++) {
	// 	double xstart = -3.0;
	// 	for (int j=0;j<n_case;j++) {
	// 		grid[i*n_case+j][0] = xstart;
	// 		grid[i*n_case+j][1] = ystart;
	// 		xstart+=(6./double(n_case));			
	// 	}
	// 	ystart+=(6./double(n_case));
	// }
	// double grain = 6.0/double(n_case);
	// for (int i=0;i<n_case*n_case;i++) {
	// 	p_goal[i] = compute_PGoal(grid[i][0], grid[i][1], 2.0, grain);
	// 	std::cout << grid[i][0] << "," << grid[i][1] << " " << p_goal[i] << std::endl;
	// }



 //  	double lower [2] = {-0.6, 1.2};
 //  	double upper [2] = {-0.4, 1.4};
 //  	double mu [2] = {-0.235, 1.193};
 //  	double cov = 0.7;
 //  	double lower2 [2];
 //  	double upper2 [2];

	// for (int i=0;i<2;i++) {
	// 	lower2[i] = (lower[i]-mu[i])/sqrt(cov);
	// 	upper2[i] = (upper[i]-mu[i])/sqrt(cov);
	// }
	// std::cout << lower2[0] << " " << lower2[1] << std::endl;
	// std::cout << upper2[0] << " " << upper2[1] << std::endl;
	// std::cout << "l " << ND2(lower[0], lower[1], 0.0) << std::endl;
	// std::cout << "u " << ND2(upper[0], upper[1], 0.0) << std::endl;
	// std::cout << ND2(upper[0], lower[1], 0.0) << std::endl;
	// std::cout << ND2(lower[0], upper[1], 0.0) << std::endl;

	// std::cout << ND2(lower[0], lower[1], 0.0)+ND2(upper[0], upper[1], 0.0)-ND2(upper[0], lower[1], 0.0)-ND2(lower[0], upper[1], 0.0) << std::endl;
	// return ND2(lower[0], lower[1], 0.0)+ND2(upper[1], upper[1], 0.0)-ND2(upper[0], lower[1], 0.0)-ND2(lower[0], upper[1], 0.0);
	// std::cout << ND2(upper2[0], upper2[1], 0.0)+ND2(lower2[0],lower2[1],0.0)-ND2(upper2[0],lower2[1],0.0)-ND2(lower2[0], upper2[1],0.0) << std::endl;
  	// std::cout << cdf_multi(lower, upper, mu, 0.6999999999) << std::endl;
	// double a, b, c, d;
	// double lower [2] = {0.1, 0.1};
	// double upper [2] = {0.3, 0.3};
	// a = ND2(lower[0], lower[1], 0.0);
	// b = ND2(lower[0], upper[1], 0.0);
	// c = ND2(upper[0], lower[1], 0.0);
	// d = ND2(upper[1], upper[1], 0.0);	
	// std::cout << a << " " << b << " " << c << " " << d <<std::endl;
	// std::cout << a+d-c-b << std::endl;
   	return 0;
}