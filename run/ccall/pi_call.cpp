#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <math.h>

using namespace std;

int dim = 2;	
int infin [2] = {2,2};
double correl = 0.0;
int maxpts = 5000;
double abseps = 0.00005;
double releps = 0;
double coeff14 = -0.5096;
double coeff31 = 0.5096;
double coeff912 = 1.0;
double coeff124 = -0.3267;
double coeff49 = -6.2527;
double coeff116 = -1.0;
double coeff63 = 6.2527;
double coeff311 = 0.3267;
double coeff3 = 1.3768;
double coeff4 = -1.3768;
double coeff6 = -3.1111;
double coeff9 = 3.1111;
double coeff11 = -0.3273;
double coeff12 = 0.3273;

extern "C" {
	int mvndst_(int *, double (*)[2], double (*)[2], int (*)[2], double *, int *,
                double *, double *, double *, double *, int *);
}
void softmax(double *p, double *v, double b) {
	double sum = 0.0;
	double tmp[5];
	// for (int i=0;i<5;i++) std::cout << v[i] << " "; std::cout << std::endl;
	for (int i=0;i<5;i++) {
		tmp[i] = exp(v[i]*b);
		sum+=tmp[i];		
	}		
	// for (int i=0;i<5;i++) std::cout << tmp[i] << " "; std::cout << std::endl;
	for (int i=0;i<5;i++) {
		p[i] = tmp[i]/sum;		
	}
	// for (int i=0;i<5;i++) std::cout << p[i] << " "; std::cout << std::endl;
	// int ind=-1;
	// for (int i=0;i<5;i++) {
	// 	if (isnan(p[i])) {
	// 		ind = i;
	// 		break;
	// 	}
	// }
	// // std::cout << ind << std::endl;
	// if (ind!=-1) {
	// 	for (int i=0;i<5;i++) {
	// 		p[i] = 0.0001;
	// 	}
	// }
	// p[ind] = 0.9996;
	// for (int i=0;i<5;i++) std::cout << p[i] << " "; std::cout << std::endl;
	// std::cout << ind << std::endl;
	// std::cout << std::endl;
	for (int i=0;i<5;i++) {
		if (p[i] == 0) {
			sum = 0.0;
			for (int i=0;i<5;i++) {
				p[i]+=1e-4;
				sum+=p[i];
			}
			for (int i=0;i<5;i++) {
				p[i]/=sum;
			}
			return;
		}
	}	
}
double sum_prod(double *a, double *b, int n) {
	double tmp = 0.0;
	for (int i=0;i<n;i++) {
		tmp+=(a[i]*b[i]);
	}
	return tmp;
}
void get_exactPosition(double *xy, int pos) {
	if (pos==0) {xy[0]=0.0;xy[1]=0.0;}
	else if (pos==1) {xy[0]=0.0;xy[1]=0.235;}
	else if (pos==2) {xy[0]=0.0;xy[1]=0.47;}
	else if (pos==3) {xy[0]=-0.19;xy[1]=0.608;}
	else if (pos==4) {xy[0]=0.19;xy[1]=0.608;}
	else if (pos==5) {xy[0]=-0.827;xy[1]=0.601;}		
	else if (pos==6) {xy[0]=-0.604;xy[1]=0.674;}						 
	else if (pos==7) {xy[0]=-0.38;xy[1]=0.746;}
	else if (pos==8) {xy[0]=0.38;xy[1]=0.746;}						 
	else if (pos==9) {xy[0]=0.604;xy[1]=0.674;}						 
	else if (pos==10) {xy[0]=0.827;xy[1]=0.601;}						 
	else if (pos==11) {xy[0]=-0.308;xy[1]=0.97;}
	else if (pos==12) {xy[0]=0.308;xy[1]=0.97;}
	else if (pos==13) {xy[0]=-0.235;xy[1]=1.193;}
	else if (pos==14) {xy[0]=0.235;xy[1]=1.193;}						 
}
double compute_Ppos(double x, double y, int pos, double varPos) {
	double meanxy [2];
	double stdev = sqrt(varPos);
	get_exactPosition(meanxy, pos);	
	double lower [2] = {(x-meanxy[0])/stdev, (y-meanxy[1])/stdev};
	double upper [2] = {(x+0.2-meanxy[0])/stdev, (y+0.2-meanxy[1])/stdev};
	// int dim = 2;	
	// int infin [2] = {2,2};
	// double correl = 0.0;
	// int maxpts = 5000;
	// double abseps = 0.00005;
	// double releps = 0;
	double error;
	double value;
	int inform;
	mvndst_(&dim, &lower, &upper, &infin, &correl, &maxpts, &abseps, &releps, &error, &value, &inform);
	return value;	
}
double compute_PGoal(double x, double y, double varGoal) {
	double meanxy [2] = {-0.235,1.193};
	double stdev = sqrt(varGoal);	
	double lower [2] = {(x-meanxy[0])/stdev, (y-meanxy[1])/stdev};
	double upper [2] = {(x+0.2-meanxy[0])/stdev, (y+0.2-meanxy[1])/stdev};
	// int dim = 2;	
	// int infin [2] = {2,2};
	// double correl = 0.0;
	// int maxpts = 5000;
	// double abseps = 0.00005;
	// double releps = 0;
	double error;
	double value;
	int inform;
	mvndst_(&dim, &lower, &upper, &infin, &correl, &maxpts, &abseps, &releps, &error, &value, &inform);
	return value;	
}
double in_action(double xp, double yp, double xt, double yt, int a, int p, int pp) {
	// std::cout << "P(x,y)="<<xp<<","<<yp << " P(test)="<<xt<<","<<yt << " action="<< a << " pos="<< p << " prev_pos=" << pp << std::endl;
	if (p==1) { // 1b
		double delta=yt-yp;
		if (pp==0 && a==0 && delta>0) return 1;
		else if (pp==2 && a==0 && delta<0) return 1;
		else if (pp==0 && a==2 && delta<0) return 1;			
		else if (pp==2 && a==2 && delta>0) return 1;
	}
	else if (p==2) {
		double d14 = yt-(xt*coeff14+(yp-coeff14*xp));
		double d31 = yt-(xt*coeff31+(yp-coeff31*xp));
		if (pp==1) {
			if (a==3 && xt > xp && d14>0) return 1;			
			else if (a==1 && xt < xp && d31>0) return 1;
			else if (a==2 && d14<0 && d31<0) return 1;
		}
		else if (pp==3) {
			if (a==1 && xt > xp && d14>0) return 1;			
			else if (a==2 && xt < xp && d31>0) return 1;
			else if (a==3 && d14<0 && d31<0) return 1;
		}
		else if (pp==4) {
			if (a==2 && xt > xp && d14>0) return 1;			
			else if (a==3 && xt < xp && d31>0) return 1;
			else if (a==1 && d14<0 && d31<0) return 1;
		}
	}
	else if (p==8) {
		double d912 = yt-(xt*coeff912+(yp-coeff912*xp));
		double d124 = yt-(xt*coeff124+(yp-coeff124*xp));
		double d49 = yt-(xt*coeff49+(yp-coeff49*xp));
		if (d49>0 && d912<0) { // first cadran
			if (pp==12 && a==1) return 1;
			else if (pp==4 && a==3) return 1;			
			else if (pp==9 && a==2) return 1;
		}
		else if (d124>0 && d912>0) { // second cadran
			if (pp==12 && a==2) return 1;
			else if (pp==4 && a==1) return 1;			
			else if (pp==9 && a==3) return 1;
		}
		else if (d124<0 && d49>0) { // third cadran
			if (pp==12 && a==3) return 1;
			else if (pp==4 && a==0) return 1;			
			else if (pp==9 && a==1) return 1;
		}
	}
	else if (p==7) {
		double d311 = yt-(xt*coeff311+(yp-coeff311*xp));
		double d116 = yt-(xt*coeff116+(yp-coeff116*xp));
		double d63 = yt-(xt*coeff63+(yp-coeff63*xp));
		if (d63<0 && d311<0) { // first cadran
			if (pp==3 && a==2) return 1;
			else if (pp==11 && a==1) return 1;			
			else if (pp==6 && a==3) return 1;
		}
		else if (d311>0 && d116>0) { // second cadran
			if (pp==3 && a==3) return 1;
			else if (pp==11 && a==2) return 1;			
			else if (pp==6 && a==1) return 1;
		}
		else if (d116<0 && d63>0) { // third cadran
			if (pp==3 && a==1) return 1;
			else if (pp==11 && a==3) return 1;			
			else if (pp==6 && a==2) return 1;
		}
	}
	else if (p==3) {
		double d = yt-(xt*coeff3+(yp-coeff3*xp));
		if (d>0) {
			if (pp==2 && a==0) return 1;
			else if (pp=7 && a==2) return 1;			
		}
		else if (d<0) {
			if (pp==2 && a==2) return 1;
			else if (pp=7 && a==0) return 1;			
		}
	}
	return 0;
}
double sferes_call(double * fit, const int N, const char* data_dir, double beta_, double gamma_, double eta_)
{	
	///////////////////
	// parameters
	double beta=0.0+beta_*(200.0-0.0);
	double gamma=0.0+gamma_*(1.0-0.0);	
	double eta=0.0+(1.0-0.0)*eta_;		

	int n_state = 3;
	int n_action = 4;
	int size_trials [N];
	int nb_points = 0;
	double varPos = gamma;
	double varGoal = 0.0;
	int n_case = 30;	
	double grid [n_case*n_case] [2];
	double p_goal [n_case*n_case];
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
	///////////////////
	const char* _data_dir = data_dir;	
	// LOADING INFO.txt
	std::string fileinfo = _data_dir;
	std::string filedata = _data_dir;
	fileinfo.append("info.txt");
	filedata.append("possarpossible.txt");
	std::ifstream data_fileinfo(fileinfo.c_str());
	string line;
	if (data_fileinfo.is_open()) {
		getline (data_fileinfo, line);
		stringstream stream(line);
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
			stringstream stream(line);
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
			xstart+=0.2;
			p_goal[i*n_case+j] = 1.0;
		}
		ystart+=0.2;
	}

	

	int index = 0;
	for (int tr=0;tr<N;tr++) {
	// for (int tr=0;tr<1;tr++) {
		// START TRIALS
		varPos = gamma;
		previous_pos = 1;
		for (int st=0;st<size_trials[tr]-1;st++) {		
			nb_possible = 0;
			int t = 0;
			for (int j=4;j<8;j++) {
				nb_possible+=sarp[index][j];
				if (sarp[index][j]==1) {
					ind_action[t] = j-4;
					t += 1;
				}			
			}		
			pos = sarp[index][0];
			state = sarp[index][1];
			action = sarp[index][2];
			reward = sarp[index][3];
			// std::cout << "trial=" << tr <<","<<st<< " pos=" << pos << " state=" << state << " action=" << action << " reward=" << reward << " nb_possible=" << nb_possible << " vPos=" << varPos << " vG=" << varGoal << std::endl;
			// COMPUTE VALUE
			if (nb_possible>1) {
				double q_values [nb_possible];
				for (int j=0;j<nb_possible;j++) {
					q_values[j] = 0.0;										
					// on regarde chaque position pour decider
					for (int k=0;k<n_case*n_case;k++){
						double subspace = 0.0;
						double Ppos = compute_Ppos(grid[k][0], grid[k][1], pos, varPos);
						// on regarde chaque position qui sont du coté de l'action
						for (int l=0;l<n_case*n_case;l++){
							double in = in_action(grid[k][0], grid[k][1], grid[l][0], grid[l][1], ind_action[j], pos, previous_pos);
							subspace+=(p_goal[l]*in);															
						}
						q_values[j] += (subspace*Ppos);
					}
					
				}			
			}			
			// UPDATE VALUE
			varPos+=gamma;
			previous_pos = pos;
			if (reward==1) {
				varGoal = (1.0-eta)*varGoal + eta * varPos;
				// FILL PGOAL
				for (int i=0;i<n_case*n_case;i++) {
					p_goal[i] = compute_PGoal(grid[i][0], grid[i][1], varGoal);
				}
			}
			index+=1;
		}
		//GUIDAGE
		if (sarp[index-1][3] == 0) {
			// std::cout << "guidage" << std::endl;
			varPos = 6.0*gamma;
			varGoal = (1.0-eta)*varGoal + eta * varPos;
		}
		index+=1;
		// std::cout << std::endl;
	}


	
	// if (isnan(fit[0]) || isinf(fit[0]) || isinf(fit[1]) || isnan(fit[1]) || fit[0]<-10000 || fit[1]<-10000) {
	// 	fit[0]=-1000.0;
	// 	fit[1]=-1000.0;
	// 	return;
	// }
	// else {
	// 	fit[0]+=2000.0;
	// 	fit[1]+=500.0;
	// 	return ;
	// }	
}
