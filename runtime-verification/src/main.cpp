#define _USE_MATH_DEFINES
#define NO_ALLOCA //The KDTree library might use the alloca function, which is not portable, nor good.
#include <cmath>
#include <iostream>
#include <ctime>
#include "Monitor.h"
#include "System.h"
#include "Property.h"
#include "RRT.h"
#include "TimedRRT.h"
#include "Plotter.h"
#include "OdeSolver.h"
#include "transient.h"
#include "hspiceInterface.h"
#include "pll.h"
#include "Monitor.h"
#include "Frequency.h"

#include "kdtree.h"
using namespace std;

#define NEW_RRT_TDO		0
#define LOAD_RRT 		1
#define NEW_RRT_PLL		2

Monitor* createMonitor_1(RRT* rrt){
	Monitor* m_const5mv = new Monitor(rrt, CONSTANT);
	m_const5mv->push_back(0.005);

	return m_const5mv;
}

//PLL locking property
// pll is locked when (v(pll_e)-v(pll_eb)) is settled to near zero for some time. 
Monitor* pllLockingPropertyMonitor(RRT* rrt){
	//(v(pll_e)-v(pll_eb))
	Monitor* m0 = new Monitor(rrt, ANALOG_BINARY_SUB );		
	m0->push_back(pll_e);
	m0->push_back(pll_eb);		

	//Norm( v(pll_e)-v(pll_eb) )
	Monitor* m1 = new Monitor( rrt, ANALOG_NORM_L2 ) ;
	m1->push_back(m0);

	//0.005mv
	Monitor* m_const5mv = new Monitor(rrt, CONSTANT);
	m_const5mv->push_back(0.005);
	
	Monitor* m3 = new Monitor(rrt, ANALOG_DIFF_SIBLING);
	m3->push_back(m1);
	
	//Norm( v(pll_e)-v(pll_eb) )<0.005mv
	Monitor* m2 = new Monitor(rrt, LOGIC_LESS_THAN_OR_EQUAL ) ;
	m2->push_back(m3);
	m2->push_back(m_const5mv);


	//Monitor* m2 = new Monitor(rrt, ANALOG_SHIFT ) ;			// m2 = x[t-(10e-7)]
	//m2->push_back(10e-7);

	//Monitor* m3 = new Monitor(rrt, ANALOG_BINARY_MUL) ;		// m3 = (x_2 + 0.5) * x[t-(10e-7)]
	//m3->push_back(m1); 
	//m3->push_back(m2);
	
	return m2 ;
}

void kernel_RRT_TDO(bool generatePlot, string outputFileName, Plotter* plotter, int iterations, double simulationTime){
	double* initialState = new double[3];

	//For oscillation:
	initialState[0] = 0.8;			//This initial condition, under low variation parameters are indicator of 
	initialState[1] = 0.04e-3;		//correct circuit that will start a healty oscillation
	initialState[2] = 0;

	//For No-Oscillation
	//initialState[0] = 0.131;			//This initial condition cause the circuit to do not oscillate. 
	//initialState[1] = 0.055;
	//initialState[2] = 0;
	// Constructing System S
	System* circuit = new System(TDO);

	TimedRRT rrt = TimedRRT(	2, //dimension
		iterations, //k
		simulationTime, //Simulation Time
		"Tunnel Diode Oscillator");
	Monitor* monitor = createMonitor_1(&rrt);
	rrt.addMonitor(monitor);
	rrt.setBound(0, -0.2, 1.2 );	//First Dimension = VC
	rrt.setBound(1, -0.02, 0.08 );  //Second Dimension = IL
	rrt.setdt(10e-9);
	rrt.setSystem(circuit);
	rrt.build(initialState, 0.005 /*variation*/);
	rrt.save(outputFileName);
	if(generatePlot) plotter->plotRRT("lines", rrt.getName().c_str(), "test", rrt, "v_C", "i_L", "t");
}

void kernel_RRT_PLL(bool generatePlot, string outputFileName, Plotter* plotter, int iterations, double simulationTime, double dt){
	int dim = 16;
	double* initialState = new double[dim+1];
	setInitialPLLState(initialState);
	// Constructing System S
	System* circuit = new System(PLL);
	TimedRRT rrt = TimedRRT(	dim, //dimension
		iterations, //k
		simulationTime, //Simulation Time
		"Tunnel Diode Oscillator");
	
	Monitor* pll_lock_property = pllLockingPropertyMonitor(&rrt);
	rrt.addMonitor(pll_lock_property);

	for(int i=0;i<dim-1;i++){
		rrt.setBound(i, -1, +1);
	}
	rrt.setdt(dt);
	rrt.setSystem(circuit);
	cout << "RRT Init." << endl ;
	rrt.build(initialState, 0.001 /*variation*/);
	cout << "RRT Constructed" << endl ;
	rrt.save(outputFileName);
	cout << "RRT Saved" << endl ;
	if(generatePlot)  plotter->plotTrace(rrt, pll_e, pll_eb, pll_time, simulationTime, dt);
}



void kernel_RRT(int mode, bool generatePlot,string inputFileName, string outputFileName, Plotter* plotter){
	if(mode==NEW_RRT_TDO){
		kernel_RRT_TDO(generatePlot, outputFileName, plotter, 100, 7e-6);
	}else if
		(mode == NEW_RRT_PLL){
		//kernel_RRT_PLL(generatePlot, outputFileName, plotter, 10, 100e-6);
		kernel_RRT_PLL(generatePlot, outputFileName, plotter, 10, 100e-6, 0.1e-6);
	}else if(mode == LOAD_RRT){
		TimedRRT rrt = TimedRRT(inputFileName);
		if(generatePlot) plotter->plotRRT("lines", rrt.getName().c_str(), "test", rrt,  "v_C", "i_L", "t");
	}
}


//TODO: 
//This function must be replaced with an actual simulation
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


<<<<<<< HEAD

void date_2013_experiments(){
=======
int main (int argc, const char * argv[]){
	srand((unsigned int)time(0));
	//Plotter* plotter = new Plotter("/Applications/Gnuplot.app/Contents/Resources/bin");
	//Plotter* plotter = new Plotter("C:\\Progra~1\\gnuplot\\bin\\gnuplot.exe -persist");
>>>>>>> origin/master
	Plotter* plotter = new Plotter("C:\\opt\\gnuplot\\bin\\gnuplot.exe -persist");

	int mode = NEW_RRT_PLL ; // NEW_RRT_TDO // LOAD_RRT // NEW_RRT_PLL
	bool generatePlot = true ;//true ;
	string inputFileName = "test.rrt";
	string outputFileName = "tdo-3000.rrt" ;

	kernel_MC();
	kernel_RRT(mode, generatePlot, inputFileName, outputFileName, plotter);
}

void fft_experiments(){
	double f0 = 1; //initial freqency
	double f1 = 10; //final freq
	double t = 5; //time interval
	int size = 1000; //number of samples
	Frequency f;
	
	vector<double> signal = f.generatoreSweepWaveform(f0, f1, t, size); //generate a sin wave sweeping from f0 to f1
	f.init();
	for(int i=0;i<size;i++){
		f.sdft();
		//if(++idx==N) idx=0;
	}
	
	
	//double powr1[N/2];
    //f.powr_spectrum(powr1);

	Plotter* plotter = new Plotter("C:\\opt\\gnuplot\\bin\\gnuplot.exe -persist");
	plotter->drawTrace(signal, t);
	plotter->saveToPdf("test.pdf");
}

void kdtree_experiment(){
	struct kdtree *kd;
    struct kdres  *set;
	
	int d=3;

	double* min = new double[d];
    double* max = new double[d];
    //default values for min-max
    for(int i=0;i<d;i++){
        min[i]=0;
        max[i]=1;
    }

	
	 kd = kd_create(d);

	 for(int i=0;i<100;i++){
		 node* q_sample = new node(d);
	     q_sample->randomize(min, max);
		 
		 kd_insert(kd,q_sample->get(), q_sample);
	 }

	node* q_sample = new node(d);
	q_sample->randomize(min, max);

	set = kd_nearest(kd, q_sample->get());
    if(kd_res_size(set)>0){
        node* res =  (node*) kd_res_item_data (set);
		cout << res->toString() << endl ;
    }else{
        cout << "[error]" << endl ;
    }

}



int main (int argc, const char * argv[]){
	srand((unsigned int)time(0));
	//fft_experiments();
	kdtree_experiment();
	system("PAUSE");
	return 0;
}


void setInitialPLLState(double* state){
	/*
	*
	* .nodeset
	+ e =  -3.5000
	+ eb =  -3.3000
	+ in =   0.
	+ inb =   0.
	+ mout =  -3.5000
	+ moutb =  -3.5000
	+ osc =-999.7550m
	+ oscb =  -1.0002
	+ out =  -3.5000
	+ outb =  -3.3000
	+ xvco.c =   1.0000
	+ xvco.s =  61.2612u
	+ xvco.s_clip =  81.6816u
	+ xpd.clip1 =   0.
	+ xpd.clip2 = 653.4528u
	+ xpd.n1 =   0.
	*/

	state[pll_e] = -3.5000 ;
	state[pll_eb] = -3.3000 ;
	state[pll_in] = 0 ;
	state[pll_inb] = 0 ;
	state[pll_mout] = -3.5000 ;
	state[pll_moutb] = -3.5000 ;
	state[pll_osc] = -999.7550e-3 ;
	state[pll_oscb] = -1.0002 ;
	state[pll_out] = -3.5000 ;
	state[pll_outb] = -3.3000 ;
	state[pll_xvco_c] = 1.0000 ;
	state[pll_xvco_s] = 61.2612e-6 ;
	state[pll_xvco_s_clip] = 81.6816e-6 ;
	state[pll_xpd_clip1] = 0 ;
	state[pll_xpd_clip2] = 653.4528e-6 ;
	state[pll_xpd_n1] = 0 ;
	state[pll_time] = 0 ;
}

