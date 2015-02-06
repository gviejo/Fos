#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iterator>
#include <math.h>
#include "normal.cpp"

// double max(double a, double b) {
// 	if (a>b)
//         return a;
// 	else
//         return b;
// }
// double min(double a, double b) {
// 	if (a>b)
//         return b;
// 	else
//         return a;
// }
// double sgn(double a) {
// 	if (a>0)
//         return 1.;
// 	else if (a<0)
//         return -1.;
// 	else
//         return 0.;
// }
// double ndf(double t) {
// 	return 0.398942280401433*exp(-t*t/2);
// }
// double nc(double x) {
// 	double result;
// 	if (x<-7.)
//         result = ndf(x)/sqrt(1.+x*x);
// 	else if (x>7.)
//         result = 1. - nc(-x);
// 	else {
// 		result = 0.2316419;
// 		double a[5] = {0.31938153,-0.356563782,1.781477937,-1.821255978,1.330274429};
// 		result=1./(1+result*fabs(x));
// 		result=1-ndf(x)*(result*(a[0]+result*(a[1]+result*(a[2]+result*(a[3]+result*a[4])))));
// 		if (x<=0.) result=1.-result;
// 	}
// 	return result;
// }
// double fxy(double x, double y, double a, double b, double rho) {
//     double a_s;
//     double b_s;
// 	double result;
//     a_s = a / sqrt(2 * (1 - rho * rho));
//     b_s = b / sqrt(2 * (1 - rho * rho));
//     result = exp(a_s * (2 * x - a_s) + b_s * (2 * y - b_s) + 2 * rho * (x - a_s) * (y - b_s));
// 	return result;
// }
// double Ntwo(double a, double b, double rho) {
//     double aij[4]={0.325303, 0.4211071, 0.1334425, 0.006374323};
//     double bij[4]={0.1337764, 0.6243247, 1.3425378,2.2626645};
//     int i;
//     int j;
//     double result;
//     result = 0;
//     for(i=0;i<=3;i++) {
//         for(j=0;j<=3;j++) {
//           	result+=aij[i] * aij[j] * fxy(bij[i], bij[j], a, b, rho); 
//         }
//     }
//     result = result * sqrt(1 - rho * rho) / M_PI;
// return result;
// }
// double ND2(double a, double b, double rho)
// {
//     double rho1;
//     double rho2;
//     double denominator;
//     double result;

//     if (a * b * rho <= 0) {
//         if (a <= 0 && b <= 0 && rho <= 0)
//             result = Ntwo(a, b, rho);
//         else if (a <= 0 && b * rho >= 0)
//             result = nc(a) - Ntwo(a, -b, -rho);
//         else if (b <= 0 && rho >= 0)
//             result = nc(b) - Ntwo(-a, b, -rho);
//         else
//             result = nc(a) + nc(b) - 1 + Ntwo(-a, -b, rho);
//     }
//     else
//             {
//         denominator = sqrt(a * a - 2 * rho * a * b + b * b);
//         rho1 = (rho * a - b) * sgn(a) / denominator;
//         rho2 = (rho * b - a) * sgn(b) / denominator;
//         result = ND2(a, 0, rho1) + ND2(b, 0, rho2) - (1 - sgn(a) * sgn(b)) / 4;
//     }
//     if (result < 0) result = 0;
//     return result;
// }


// double cdf_multi(double *lower, double *upper, double *mu, double cov) {
// 	for (int i=0;i<2;i++) {
// 		lower[i] = (lower[i]-mu[i])/sqrt(cov);
// 		upper[i] = (upper[i]-mu[i])/sqrt(cov);
// 	}
// 	return bivnor(lower[0], lower[1], 0.0)+bivnor(upper[1], upper[1], 0.0)-bivnor(upper[0], lower[1], 0.0)-bivnor(lower[0], upper[1], 0.0);

// }
void softmax(double *p, double *v, double beta, int n) {
	double sum = 0.0;
	double tmp[n];
	int have_inf=0;
	for (int i=0;i<n;i++) {
		tmp[i] = exp(v[i]*beta);
		if (isinf(tmp[i])) {
			have_inf = 1;
			tmp[i]=1.0-double(n)*1e-8;
		}
		sum+=tmp[i];		
	}			
	if (have_inf==1) {
		for (int i=0;i<n;i++) tmp[i] += 1e-8;			
	}
	int have_zero=0;
	for (int i=0;i<n;i++) {
		p[i] = tmp[i]/sum;
		if (p[i]==0.0) have_zero=1;
	}
	if (have_zero==1) {
		sum = 0.0;
		for (int i=0;i<n;i++) {
			p[i] += 1e-8;
			sum+=p[i];
		}
		for (int i=0;i<n;i++) {
			p[i] = p[i]/sum;			
		}
	}
}
void get_exactPosition(double *xy, int pos) {
	// if (pos==0) {xy[0]=0.0;xy[1]=0.0;}
	if (pos==1) {xy[0]=0.0;xy[1]=0.235;}
	else if (pos==2) {xy[0]=0.0;xy[1]=0.47;}
	else if (pos==3) {xy[0]=-0.19;xy[1]=0.608;}
	else if (pos==4) {xy[0]=0.19;xy[1]=0.608;}
	// else if (pos==5) {xy[0]=-0.827;xy[1]=0.601;}		
	else if (pos==6) {xy[0]=-0.604;xy[1]=0.674;}						 
	else if (pos==7) {xy[0]=-0.38;xy[1]=0.746;}
	else if (pos==8) {xy[0]=0.38;xy[1]=0.746;}						 
	else if (pos==9) {xy[0]=0.604;xy[1]=0.674;}						 
	// else if (pos==10) {xy[0]=0.827;xy[1]=0.601;}						 
	else if (pos==11) {xy[0]=-0.308;xy[1]=0.97;}
	else if (pos==12) {xy[0]=0.308;xy[1]=0.97;}
	// else if (pos==13) {xy[0]=-0.235;xy[1]=1.193;}
	// else if (pos==14) {xy[0]=0.235;xy[1]=1.193;}						 
}
double compute_Ppos(double x, double y, int pos, double varPos, double grain_) {	
	double meanxy [2];
	get_exactPosition(meanxy, pos);
	double stdev = sqrt(varPos);			
	double lower [2] = {(x-meanxy[0])/stdev, (y-meanxy[1])/stdev};
	double upper [2] = {(x+grain_-meanxy[0])/stdev, (y+grain_-meanxy[1])/stdev};
	return bivnor(upper[0], upper[1], 0.0)+bivnor(lower[0],lower[1],0.0)-bivnor(upper[0],lower[1],0.0)-bivnor(lower[0], upper[1],0.0);
}
double compute_PGoal(double x, double y, double varGoal, double grain_) {
	double stdev = sqrt(varGoal);
	double lower [2] = {(x-(-0.235))/stdev, (y-1.193)/stdev};
	double upper [2] = {(x+grain_-(-0.235))/stdev, (y+grain_-1.193)/stdev};
	return bivnor(upper[0], upper[1], 0.0)+bivnor(lower[0],lower[1],0.0)-bivnor(upper[0],lower[1],0.0)-bivnor(lower[0], upper[1],0.0);
}
int in_Yaction(double xp, double yp, double xt, double yt, int p) {
	double coeff14 = -0.5096;
	double coeff31 = 0.5096;
	double coeff912 = 1.0;
	double coeff124 = -0.3267;
	double coeff49 = -6.2527;
	double coeff116 = -1.0;
	double coeff63 = 6.2527;
	double coeff311 = 0.3267;
	if (p==2) {
		double d14 = yt-(xt*coeff14+(yp-coeff14*xp));
		double d31 = yt-(xt*coeff31+(yp-coeff31*xp));
		if (xt > xp && d14>=0) {
			return 0;
		}
		else if (xt < xp && d31>=0) {
			return 1;
		}
		else if (d14<0 && d31<0) {
			return 2;
		}
		else return -1;
	}
	else if (p==8) {
		double d912 = yt-(xt*coeff912+(yp-coeff912*xp));
		double d124 = yt-(xt*coeff124+(yp-coeff124*xp));
		double d49 = yt-(xt*coeff49+(yp-coeff49*xp));
		// std::cout << d912 << " " << d124 << " " << d49 << std::endl;
		if (d49>=0 && d912<0) { // first cadran
			return 0;
		}
		else if (d124>=0 && d912>=0) { // second cadran			
			return 1;
		}
		else if (d124<0 && d49<0) { // third cadran
			return 2;
		}
		else return -1;
	}
	else if (p==7) {
		double d311 = yt-(xt*coeff311+(yp-coeff311*xp));
		double d116 = yt-(xt*coeff116+(yp-coeff116*xp));
		double d63 = yt-(xt*coeff63+(yp-coeff63*xp));
		if (d63<0 && d311<0) { // first cadran
			return 0;
		}
		else if (d311>=0 && d116>=0) { // second cadran
			return 1;
		}
		else if (d116<0 && d63>=0) { // third cadran
			return 2;
		}
		else return -1;
	}	
	return -1;
}
int in_Iaction(double xp, double yp, double xt, double yt, int p, double *coeff, int *ind_pgi) {
	double d = yt-(xt*coeff[ind_pgi[p]]+(yp-coeff[ind_pgi[p]]*xp));
	if (d>=0.0) return 0;
	else return 1;	
}
void update_goal(double *pgoal, double (*PG8) [3], double (*PG7) [3], double (*PG2) [3], double (*PGI) [2][7], double varGoal, double (*grid) [2], double *coeff, int *ind_pgi) {	
	double grain_ = 6.0/30.0;
	// double tmp [31][31];
	for (int i=0;i<900;i++) {
		for (int j=0;j<3;j++) {
			PG8[i][j] = 0.0;
			PG7[i][j] = 0.0;
			PG2[i][j] = 0.0;
		}
		for (int j=0;j<7;j++) {
			PGI[i][0][j] = 0.0;
			PGI[i][1][j] = 0.0;
		}
		// prefill_pgoal[i] = compute_prefillPGoal(grid[i][0], grid[i][1], varGoal, grain_);
	}
	int ind_cadran;
	int ind_pos2 [7] = {1,3,4,6,9,11,12};	
	for (int i=0;i<900;i++) {		
		pgoal[i] = compute_PGoal(grid[i][0], grid[i][1], varGoal, grain_);
		for (int j=0;j<900;j++) {
			ind_cadran = in_Yaction(grid[j][0], grid[j][1], grid[i][0], grid[i][1], 8);
			if (ind_cadran!=-1) PG8[j][ind_cadran] += pgoal[i];
			ind_cadran = in_Yaction(grid[j][0], grid[j][1], grid[i][0], grid[i][1], 7);
			if (ind_cadran!=-1) PG7[j][ind_cadran] += pgoal[i];
			ind_cadran = in_Yaction(grid[j][0], grid[j][1], grid[i][0], grid[i][1], 2);
			if (ind_cadran!=-1) PG2[j][ind_cadran] += pgoal[i];

			for (int k=0;k<7;k++) {				
				ind_cadran = in_Iaction(grid[j][0], grid[j][1], grid[i][0], grid[i][1], ind_pos2[k], coeff, ind_pgi);
				PGI[j][ind_cadran][k] += pgoal[i];
			}
		}
	}
}
void get_Yorder_action(int *order, int p, int pp) {
	if (p==2) {
		if (pp==1) {
			order[0] = 2;
			order[1] = 0;
			order[2] = 1;
			return ;
		}
		else if (pp==3) {
			order[0] = 0;
			order[1] = 1;
			order[2] = 2;
			return ;
		}
		else if (pp==4) {
			order[0] = 1;
			order[1] = 2;
			order[2] = 0;
			return ;
		}
	}
	else if (p==8) {		
		if (pp==12) {
			order[0] = 0;
			order[1] = 1;
			order[2] = 2;
			return ;
		}
		else if (pp==4) {
			order[0] = 2;
			order[1] = 0;
			order[2] = 1;
			return ;
		}
		else if (pp==9) {
			order[0] = 1;
			order[1] = 2;
			order[2] = 0;
			return ;
		}
	}
	else if (p==7) {		
		if (pp==3) {
			order[0] = 1;
			order[1] = 2;
			order[2] = 0;
			return ;
		}
		else if (pp==11) {
			order[0] = 0;
			order[1] = 1;
			order[2] = 2;
			return ;
		}
		else if (pp==6) {
			order[0] = 2;
			order[1] = 0;
			order[2] = 1;
			return ;
		}
	}
}
void get_Iorder_action(int *order, int p, int pp) {
	if (p==1 && pp==0) {order[0]=0;order[1]=1;return;}
	else if (pp==2 && p==4) {order[0]=0;order[1]=1;return;}
	else if (pp==2 && p==3) {order[0]=0;order[1]=1;return;}
	else if (pp==5 && p==6) {order[0]=0;order[1]=1;return;}
	else if (pp==7 && p==11) {order[0]=0;order[1]=1;return;}
	else if (pp==8 && p==12) {order[0]=0;order[1]=1;return;}
	else if (pp==10 && p==9) {order[0]=0;order[1]=1;return;}
	else {		
		order[0]=1;
		order[1]=0;		
		return;}
}
void compute_3qv(double *q_values, int pos, int previous_pos, double (*PG) [3], double varPos, double (*grid)[2]) {
	double grain_ = 6.0/30.0;
	int order [3];
	get_Yorder_action(order, pos, previous_pos);	
	for (int j=0;j<900;j++) {
		for (int i=0;i<3;i++) {	
			q_values[order[i]] += (PG[j][i]*compute_Ppos(grid[j][0], grid[j][1], pos, varPos, grain_));
		}
	}
}
void compute_2qv(double *q_values, int pos, int previous_pos, double (*PG) [2][7], double varPos, double (*grid)[2], int *ind_pgi) {	
	double grain_ = 6.0/30.0;
	int order [2];	
	int ind_pos = ind_pgi[pos];
	get_Iorder_action(order, pos, previous_pos);				
	for (int j=0;j<900;j++) {
		for (int i=0;i<2;i++) {			
			q_values[order[i]] += (PG[j][i][ind_pos]*compute_Ppos(grid[j][0], grid[j][1], pos, varPos, grain_));
		}
	}
}
void sferes_call(double * fit, const int N, const char* data_dir, double beta_, double gamma_, double eta_)
{	
	///////////////////
	// parameters
	double beta=0.0+beta_*(200.0-0.0);
	double gamma=0.0+gamma_*(0.999999999-0.0);	
	double eta=0.1+(0.1-0.01)*eta_;

	// std::cout << beta << " " << gamma << " " << eta << std::endl;
	const int n_state = 3;
	const int n_action = 4;
	const int n_case = 30;
	const int n_case2 = 900;
	int size_trials [N];
	int nb_points = 0;
	double varPos = gamma;
	double varGoal = 0.0;	
	double grid [n_case2] [2];
	double p_goal [n_case2];
	int reward_position = 13;
	double logLikelihood = 0.0;	
	int possible[n_action];
	int nb_possible;
	int state;
	int action;
	int reward;
	int pos;
	int previous_pos = 0;
	int ind_action [3];
	double PG8[n_case2][3];
	double PG2[n_case2][3];
	double PG7[n_case2][3];
	double PGI[n_case2][2][7];
	double q_values[4];	
	double p_a[3];	
	double coeff [7] = {0.0, 1.3768, -1.3768, -3.1111, 3.1111, -0.3273, 0.3273};
	int ind_pgi [13] = {0,0,0,1,2,0,3,0,0,4,0,5,6};
	int ind_pos2 [7] = {1,3,4,6,9,11,12};	
	// double prefill_pgoal[n_case2];
	// double prefill_ppos[n_case2];
	for (int i=0;i<n_case2;i++) {
		for (int j=0;j<3;j++) {
			PG8[i][j] = 0.0;
			PG2[i][j] = 0.0;
			PG7[i][j] = 0.0;
		}
		for (int j=0;j<7;j++) {
			PGI[i][0][j] = 0.0;
			PGI[i][1][j] = 0.0;
		}
	}
	///////////////////
	const char* _data_dir = data_dir;	
	// LOADING INFO.txt
	std::string fileinfo = _data_dir;
	std::string filedata = _data_dir;
	fileinfo.append("info.txt");
	filedata.append("possarpossible.txt");	
	std::ifstream data_fileinfo(fileinfo.c_str());
	std::string line;
	if (data_fileinfo.is_open()) {
		getline (data_fileinfo, line);
		std::stringstream stream(line);
		std::vector<int> values(
			(std::istream_iterator<int>(stream)),
			(std::istream_iterator<int>()));
		for (int j=0;j<N;j++) {
			size_trials[j] = values[j];
			nb_points+= size_trials[j];
		}
	}
	data_fileinfo.close();	
	//////////////////	
	// LOADING possarpossible.txt
	int sarp [nb_points][8];		
	std::ifstream data_filedata(filedata.c_str());
	if (data_filedata.is_open())
	{ 
		for (int i=0;i<nb_points;i++) 
		{  			
			getline (data_filedata,line);			
			std::stringstream stream(line);
			std::vector<int> values(
     			(std::istream_iterator<int>(stream)),
     			(std::istream_iterator<int>()));
			for (int j=0;j<8;j++)
			{				
				sarp[i][j] = values[j];
			}
		}
	data_filedata.close();	
	}	
	// fill x, y grid
	double ystart = -3.0;
	for (int i=0;i<n_case;i++) {
		double xstart = -3.0;
		for (int j=0;j<n_case;j++) {
			grid[i*n_case+j][0] = xstart;
			grid[i*n_case+j][1] = ystart;
			xstart+=(6./double(n_case));
			p_goal[i*n_case+j] = 0.0;			
		}
		ystart+=(6./double(n_case));
	}

	// First trial  // NO COMPUTE VALUE
	int index = 0;
	for (int st=0;st<size_trials[0]-1;st++) {				
		pos = sarp[index][0];
		state = sarp[index][1];
		action = sarp[index][2];
		reward = sarp[index][3];			
		int t = 0;		
		nb_possible = 0;
		for (int j=4;j<8;j++) {
			nb_possible+=sarp[index][j];
			if (sarp[index][j]==1) {
				ind_action[t] = j-4;					
				t+=1;
			}			
			q_values[j-4] = 0.0;
		}
		softmax(p_a, q_values, beta, nb_possible);			
		// ADDING LOGLIKELIHOOD 
		for (int i=0;i<nb_possible;i++) {					
			if (ind_action[i] == action) {
				logLikelihood += log(p_a[i]);
				// this_log += log(p_a[i]);
			}
		}		
		// UPDATE VALUE
		varPos+=gamma;
		previous_pos = pos;			
		if (reward==1) {				
			varGoal = (1.0-eta)*varGoal + eta * varPos;				
			// FILL PGOAL and PG*				
			update_goal(p_goal, PG8, PG7, PG2, PGI, varGoal, grid, coeff, ind_pgi);
		}
		index+=1;			
	}
	//GUIDAGE
	if (sarp[index-1][3] == 0) {			
		varPos = 7.0*gamma;
		varGoal = (1.0-eta)*varGoal + eta * varPos;			
		// FILL PGOAL and PG*
		update_goal(p_goal, PG8, PG7, PG2, PGI, varGoal, grid, coeff, ind_pgi);
	}
	index+=1;
	// OTHERS TRIALS 
	for (int tr=1;tr<N;tr++) {
	// for (int tr=0;tr<1;tr++) {
		// START TRIALS
		varPos = gamma;
		previous_pos = 0;
		// double this_log = 0.0;
		for (int st=0;st<size_trials[tr]-1;st++) {				
			pos = sarp[index][0];
			state = sarp[index][1];
			action = sarp[index][2];
			reward = sarp[index][3];			
			int t = 0;		
			nb_possible = 0;
			for (int j=4;j<8;j++) {
				nb_possible+=sarp[index][j];
				if (sarp[index][j]==1) {
					ind_action[t] = j-4;					
					t+=1;
				}			
				q_values[j-4] = 0.0;
			}						
			// COMPUTE VALUE
			if (nb_possible==2) {
				compute_2qv(q_values, pos, previous_pos, PGI, varPos, grid, ind_pgi);
				softmax(p_a, q_values, beta, nb_possible);			
				// ADDING LOGLIKELIHOOD 
				for (int i=0;i<nb_possible;i++) {					
					if (ind_action[i] == action) {
						logLikelihood += log(p_a[i]);
						// this_log += log(p_a[i]);
					}
				}	
				// if (tr == 3) {
				// 	std::cout << q_values[0] << " " << q_values[1] << " " << std::endl;			
				// }			
			}
			else if (nb_possible == 3) {				
				if (pos==8) compute_3qv(q_values, 8, previous_pos, PG8, varPos, grid);
				else if (pos==7) compute_3qv(q_values, 7, previous_pos, PG7, varPos, grid);
				else if (pos==2) compute_3qv(q_values, 2, previous_pos, PG2, varPos, grid);
				softmax(p_a, q_values, beta, nb_possible);
				// ADDING LOGLIKELIHOOD 
				for (int i=0;i<nb_possible;i++) {					
					if (ind_action[i] == action) {
						logLikelihood += log(p_a[i]);
						// this_log += log(p_a[i]);
					}
				}
				// if (tr == 3) {
				// 	for (int i=0;i<3;i++) { std::cout << q_values[i] << " " ;}
				// 	std::cout << std::endl;
				// }
			}
				
			// UPDATE VALUE
			varPos+=gamma;
			previous_pos = pos;			
			if (reward==1) {				
				varGoal = (1.0-eta)*varGoal + eta * varPos;				
				// FILL PGOAL and PG*				
				update_goal(p_goal, PG8, PG7, PG2, PGI, varGoal, grid, coeff, ind_pgi);
			}
			index+=1;			
		}
		//GUIDAGE
		if (sarp[index-1][3] == 0) {			
			varPos = 7.0*gamma;
			varGoal = (1.0-eta)*varGoal + eta * varPos;			
			// FILL PGOAL and PG*
			update_goal(p_goal, PG8, PG7, PG2, PGI, varGoal, grid, coeff, ind_pgi);
		}
		index+=1;
		// std::cout << tr << " " << this_log << std::endl;
		// std::cout << this_log << std::endl;

	}	

	fit[0] = logLikelihood;	

	if (isnan(fit[0]) || isinf(fit[0]) || fit[0] == 0) {
	 	fit[0]=-10000.0;		
	 	return;
	}
}

