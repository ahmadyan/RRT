#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <ctime>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <direct.h>

#include "Monitor.h"
#include "System.h"
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
#define LOAD_RRT_PLL	3
#define LOAD_RRT_TDO	4	

void kernel_RRT_TDO(vector<Monitor*> monitors, bool generatePlot, string outputFileName, Plotter* plotter, int iterations, double simulationTime, double dt){
	double* initialState = new double[3];
	int variations = 2 ;
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

	TimedRRT rrt = TimedRRT(	
		2, //dimension
		iterations, //k
		variations,
		simulationTime, //Simulation Time
		"Tunnel Diode Oscillator");

	for(int i=0;i<(int)(monitors.size());i++){
		rrt.addMonitor(monitors[i]);
	}
	rrt.setBound(0, -0.2, 1.2 );	//First Dimension = VC
	rrt.setBound(1, -0.02, 0.08 );  //Second Dimension = IL
	rrt.setdt(dt);
	rrt.setSystem(circuit);
	rrt.build(initialState, 0.1 /*variation*/);
	//rrt.simulate(initialState, 0.005 /*variation*/);  //for just simulating the circuit
	rrt.save(outputFileName);
	if(generatePlot) plotter->plotRRT("lines", rrt.getName().c_str(), "test", rrt, "v_C", "i_L", "t");
}

void kernel_RRT_PLL(vector<Monitor*> monitors, bool generatePlot, string outputFileName, Plotter* plotter, int iterations, double simulationTime, double dt){
	int dim = 16;
	int variations = 1; 
	double* initialState = new double[dim+1];
	setInitialPLLState(initialState);
	// Constructing System S
	System* circuit = new System(PLL);
	TimedRRT rrt = TimedRRT(	dim, //dimension
		iterations, //k
		variations,
		simulationTime, //Simulation Time
		"Tunnel Diode Oscillator");
	
	for(int i=0;i<(int)(monitors.size());i++){
		rrt.addMonitor(monitors[i]);
	}
	
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
	//plotter->plotTrace(rrt, pll_e, -1, pll_time, simulationTime, dt);
}

void kernel_RRT(vector<Monitor*> monitors, int mode, bool generatePlot,string inputFileName, string outputFileName, Plotter* plotter){
	if(mode==NEW_RRT_TDO){
		kernel_RRT_TDO(monitors, generatePlot, outputFileName, plotter, 300, 7e-6, 1e-9);
	}else if(mode == NEW_RRT_PLL){
		kernel_RRT_PLL(monitors, generatePlot, outputFileName, plotter, 10, 100e-6, 0.01e-6);
	}else if(mode == LOAD_RRT_PLL){
		TimedRRT rrt = TimedRRT(inputFileName);
		if(generatePlot)  plotter->plotTrace(rrt, pll_e, pll_eb, pll_time, 100e-6, 0.01e-6);
	}else if(mode == LOAD_RRT){
		TimedRRT rrt = TimedRRT(inputFileName);
		if(generatePlot) plotter->plotRRT("lines", rrt.getName().c_str(), "test", rrt,  "v_C", "i_L", "t");
	}
}

void date_2013_experiments(){
	vector<Monitor*> monitors;	//for this experiment, the monitors are empty. 
	Plotter* plotter = new Plotter("C:\\opt\\gnuplot\\bin\\gnuplot.exe -persist");
	int mode = NEW_RRT_TDO ; // NEW_RRT_TDO // LOAD_RRT // NEW_RRT_PLL, LOAD_RRT_PLL
	bool generatePlot = true ;//true ;
	string inputFileName = "test2.rrt";
	string outputFileName = "exp4_rrt_oscillation_bad_3000.rrt" ;
	//the new kernel has the monitors argument
	kernel_RRT(monitors, mode, generatePlot, inputFileName, outputFileName, plotter);
}

//	Computing the joint time-frequency space instead of only time-augmented RRT
//	work-in-progress, todo
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

	//getting the nearest node
	//struct kdres *kd_nearest(struct kdtree *tree, const double *pos);
	set = kd_nearest(kd, q_sample->get());
    if(kd_res_size(set)>0){
        node* res =  (node*) kd_res_item_data (set);
		cout << res->toString() << endl ;
    }else{
        cout << "[error]" << endl ;
    }


	//getting the set of nearest node for different time-frames
	
	//struct kdres *kd_nearest_range(struct kdtree *tree, const double *pos, double range);
	/*
	struct kdres *presults;
	presults = kd_nearest_range(kd, q_sample->get(), 0.01); 
	if(kd_res_size(presults)>0){	
		while( !kd_res_end( presults ) ) {
		// get the data and position of the current result item 
		node* res =  (node*) kd_res_item( presults, pos );
		cout << res->toString() << endl ;
		// compute the distance of the current result from the pt 
		dist = sqrt( dist_sq( pt, pos, 3 ) );

		// print out the retrieved data 
		printf( "node at (%.3f, %.3f, %.3f) is %.3f away and has data=%c\n", 
	    pos[0], pos[1], pos[2], dist, *pch );

		// go to the next entry 
		kd_res_next( presults );
		}
    }else{
        cout << "[error] 2" << endl ;
    }*/

}

#include "ParseTree.h"
void MonitorExperiment(){
	vector<Monitor*> monitors;
	Monitor* monitor = new Monitor();
	Property* pt = new Property(1);
	monitor->setProperty(pt);
	monitors.push_back(monitor);

	Plotter* plotter = new Plotter("C:\\opt\\gnuplot\\bin\\gnuplot.exe -persist");
	int mode = NEW_RRT_PLL ; // NEW_RRT_TDO // LOAD_RRT // NEW_RRT_PLL, LOAD_RRT_PLL
	bool generatePlot = true ;//true ;
	string inputFileName = "test2.rrt";
	string outputFileName = "test.rrt" ;
	kernel_RRT(monitors, mode, generatePlot, inputFileName, outputFileName, plotter);


	//pt->parseFormula("");
	//pt->printParseTree(pt->getRoot());
	//setup an example tree execution to check the monitor
}

int main (int argc, const char * argv[]){
	srand((unsigned int)time(0));
	char full[_MAX_PATH];
	_fullpath(full, ".\\", _MAX_PATH);
	cout << "Current working directory is:" << full << endl;

	//fft_experiments();
	//kdtree_experiment();
	//MonitorExperiment();
	date_2013_experiments();
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

