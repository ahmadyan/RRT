// I don't remember why this code exists in this project, It was in the main file.
// I moved it to here for clean-up. Removing it doesn't break anything either. Was this part of generating a figure or something?
// Apparantly it computes value of pi (3.14) using MC. Why would I need to do this?
// Maybe RRT vs MC?

#include <cmath>
#include <iostream>
#include <ctime>
using namespace std;

bool MC_run_simulation(){
	double x = (double)rand()/RAND_MAX;
	double  y = (double)rand()/RAND_MAX;
	double  z = x*x+y*y;
	if (z<=1){
		return 1;
	}else{
		return 0;
	}
}

void kernel_MC(){
	double pi = 0 ;
	long count=0;
	long trials = 0 ;
	long max_iteration = 100000;
	double error_tolerance = 0.05;
	double sim_epp, sim_var;
	double z_alpha_half = 2.576;
	int block_iteration = 100;
	do{
		for(int i=0;i<block_iteration;i++){
			count += MC_run_simulation();
		}
		trials += block_iteration ;
		sim_epp  = count / (double)trials;
		sim_var = sqrt((sim_epp * (1 - sim_epp))/(double)trials);
		if (sim_var==0) sim_var = sqrt(((1.0/(double)trials) * (1 - (1.0/(double)trials)))/(double)trials);
	} while ((trials<max_iteration) && ((z_alpha_half*sim_var) > error_tolerance));
	cout << trials << " " << count << " " << sim_epp << " " << sim_var << endl ;
	pi = count / (double)trials*4 ;
	cout << "pi=" << pi << endl ;
}

