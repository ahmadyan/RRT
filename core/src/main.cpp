#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <ctime>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <direct.h>
#include <random>
#include "config.h"
#include "Monitor.h"
#include "RRT.h"
#include "TimedRRT.h"
#include "Plotter.h"
#include "OdeSolver.h"
#include "transient.h"

#include "Monitor.h"
#include "Frequency.h"
#include "kdtree.h"

#include "System.h"
#include "spice.h"
#include "vanderpol.h"
#include "pll.h"
#include "tdo.h"
#include "inverter.h"
#include "josephson.h"

using namespace std;

#define NEW_RRT_TDO		0
#define LOAD_RRT 		1
#define NEW_RRT_PLL		2
#define LOAD_RRT_PLL	3
#define LOAD_RRT_TDO	4
#define NEW_RRT_INV		5
#define LOAD_RRT_INV	6
#define SIM_TDO			7
#define SIM_PLL			8
#define SIM_INV			9
#define SIM_RING		10
/*
void kernel_RRT_TDO(vector<Monitor*> monitors, bool generatePlot, string outputFileName, Plotter* plotter, int iterations, double simulationTime, double dt, bool mc){
	double* initialState = new double[3];
	int variations = 2 ;
	
	TDO* circuit = new TDO();
	//For oscillation:
	initialState[0] = 0.8;			//This initial condition, under low variation parameters are indicator of 
	initialState[1] = 0.04e-3;		//correct circuit that will start a healty oscillation
	initialState[2] = 0;

	//For No-Oscillation
	//initialState[0] = 0.131;			//This initial condition cause the circuit to do not oscillate. 
	//initialState[1] = 0.055;
	//initialState[2] = 0;
	
	int iter=0;
	if (mc){
		iter = simulationTime / dt;
	}
	else{
		iter = iterations;
	}
	TimedRRT rrt = TimedRRT(	
		2, //dimension
		iter, //k
		variations,
		simulationTime, //Simulation Time
		"Tunnel Diode Oscillator");

	for(int i=0;i<(int)(monitors.size());i++){
		rrt.addMonitor(monitors[i]);
	}
	rrt.setBound(0, -0.2, 1.2 );	//First Dimension = VC
	rrt.setBound(1, -0.02, 0.08 );  //Second Dimension = IL

	rrt.setVariationBound(0, -0.005, 0.005);	//p0, v = 300mv +- p0
	rrt.setVariationBound(1, -0.005, 0.0005);//p1, i = id(vc) +- p1
	rrt.setdt(dt);
	
	rrt.setSystem(circuit);
	if (mc){
		rrt.simulate(initialState);  //for just simulating the circuit
	}else{
		rrt.build(initialState);
	}
	rrt.save(outputFileName);
	if(generatePlot) plotter->plotRRT("lines", rrt.getName().c_str(), "test", rrt, "v_C", "i_L", "t");
}

void kernel_RRT_PLL(vector<Monitor*> monitors, bool generatePlot, string outputFileName, Plotter* plotter, int iterations, double simulationTime, double dt, bool mc){
	int dim = 16;
	int variations = 1; 
	
	double* initialState = new double[dim+1];
	PLL* circuit = new PLL();		// Constructing System S
	circuit->setInitialPLLState(initialState);

	TimedRRT rrt = TimedRRT(	dim, //dimension
		iterations, //k
		variations,
		simulationTime, //Simulation Time
		"Phased Locked Loop");
	
	for(int i=0;i<(int)(monitors.size());i++){
		rrt.addMonitor(monitors[i]);
	}
	
	for(int i=0;i<dim-1;i++){
		rrt.setBound(i, -1, +1);
	}
	rrt.setVariationBound(0, -0.001, 0.001);
	
	rrt.setdt(dt);
	rrt.setSystem(circuit);
	
	if (mc){
		rrt.simulate(initialState);  //for just simulating the circuit
	}
	else{
		rrt.build(initialState);
	}
	rrt.save(outputFileName);
	if(generatePlot)  plotter->plotTrace(rrt, pll_e, pll_eb, pll_time, simulationTime, dt);
}


void kernel_RRT(vector<Monitor*> monitors, int mode, bool generatePlot,string inputFileName, string outputFileName, Plotter* plotter){
	if(mode==NEW_RRT_TDO){
		int iterations = 10000; 
		double simTime = 2e-6;
		double dt = 5e-9; //hard-coded in hspice code
		bool monteCarlo = false; //true for rrt, false for mc
		kernel_RRT_TDO(monitors, generatePlot, outputFileName, plotter, iterations, simTime, dt, monteCarlo);
	}else if(mode == NEW_RRT_PLL){
		int iterations = 10000;
		double simTime = 100e-6;
		double dt = 0.01e-6; //hard-coded in hspice code
		bool monteCarlo = false; //true for rrt, false for mc
		kernel_RRT_PLL(monitors, generatePlot, outputFileName, plotter, iterations, simTime, dt, monteCarlo);
	}else if(mode == LOAD_RRT_PLL){
		TimedRRT rrt = TimedRRT(inputFileName);
		if(generatePlot)  plotter->plotTrace(rrt, pll_e, pll_eb, pll_time, 100e-6, 0.01e-6);
	}else if(mode == LOAD_RRT){
		TimedRRT rrt = TimedRRT(inputFileName);
		if(generatePlot) plotter->plotRRT("lines", rrt.getName().c_str(), "test", rrt,  "v_C", "i_L", "t");
	}else if (mode == NEW_RRT_INV){
		kernel_RRT_INV(monitors, generatePlot, outputFileName, plotter, 100, 1e-9, 5e-12, false);
	}else if (mode == SIM_TDO){
		int iterations = -1;
		double simTime = 5e-6;
		double dt = 5e-9; //hard-coded in hspice code
		bool monteCarlo = true; //true for rrt, false for mc
		kernel_RRT_TDO(monitors, generatePlot, outputFileName, plotter, iterations, simTime, dt, monteCarlo);
	}else if (mode == SIM_PLL){
		double simTime = 1e-4;
		double dt = 1e-7;
		int iterations = simTime / dt;
		bool monteCarlo = true; //true for rrt, false for mc
		kernel_RRT_PLL(monitors, generatePlot, outputFileName, plotter, iterations, simTime, dt, monteCarlo);
	}
}
*/

System* systemSelector(Configuration* config){
	System* circuit;
	if (config->checkParameter("edu.uiuc.csl.system.name", "Vanderpol")){
		circuit = new Vanderpol(config);
	}else if (config->checkParameter("edu.uiuc.csl.system.name", "Josephson")){
		circuit = new Josephson(config);
	}else if (config->checkParameter("edu.uiuc.csl.system.name", "Inverter")){
		circuit = new Inverter(config);
	}else if (config->checkParameter("edu.uiuc.csl.system.name", "TDO")){
		circuit = new TDO(config);
	}else if (config->checkParameter("edu.uiuc.csl.system.name", "PLL")){
		circuit = new PLL(config);
	}else if (config->checkParameter("edu.uiuc.csl.system.name", "Josephson")){
		circuit = new Josephson(config);
	}else if (config->checkParameter("edu.uiuc.csl.system.name", "Josephson")){
		circuit = new Josephson(config);
	}else{
		if (config->checkParameter("edu.uiuc.csl.system.simulator", "hspice")){
			circuit = new SPICE(config);
		}else{
			string name; config->getParameter("edu.uiuc.csl.system.name", &name);
			cout << "System " << name << " is not supported. " << endl;
		}
		
	}
	return circuit;
}

void kernel_RRT(Configuration* config){
	string outputFileName;  config->getParameter("edu.uiuc.crhc.core.options.outputFileName", &outputFileName);
	string inputFileName;  config->getParameter("edu.uiuc.crhc.core.options.inputFileName", &inputFileName);

	double simulationTime; config->getParameter("edu.uiuc.csl.core.simulation.simulationTime", &simulationTime);
	double dt; config->getParameter("edu.uiuc.csl.core.simulation.dt", &dt);
	string name; config->getParameter("edu.uiuc.csl.system.name", &name);

	int dim; config->getParameter("edu.uiuc.csl.system.dimension", &dim);
	int param; config->getParameter("edu.uiuc.csl.system.parameters", &param);
	
	System* circuit = systemSelector(config);
	
	double* initialState = new double[dim + 1]; // counting the time-augmentation as another dimension
	initialState[dim] = 0;	//time always begins at 0
	if (config->checkParameter("edu.uiuc.csl.system.ic", "file")){
		double* tmp = circuit->parseICFile(config->get("edu.uiuc.csl.system.ic.file"));
		for (int i = 0; i < dim; i++){
			initialState[i] = tmp[i];
		}
		delete tmp;
	}
	else{
		for (int i = 0; i < dim; i++){
			double ic; config->getParameter("edu.uiuc.csl.system.var.ic", i, &ic);
			initialState[i] = ic;
		}
	}
	
	int iter = 0;
	if (config->checkParameter("edu.uiuc.crhc.core.mode", "mc")){
		iter = simulationTime / dt;
	}else{
		config->getParameter("edu.uiuc.csl.core.sampling.iteration", &iter);
	}

	//todo: add this throught logbook
	cout << "----------------------------------------------" << endl;
	cout << "RRT Options: " << endl;
	cout << "Iter=" << iter << endl;
	cout << "sim-time=" << simulationTime << endl;
	cout << "dt=" << dt << endl;

	cout << "----------------------------------------------" << endl;
	TimedRRT rrt = TimedRRT(config, dim, iter, param, simulationTime, name);

	//for (int i = 0; i<(int)(monitors.size()); i++){
	//	rrt.addMonitor(monitors[i]);
	//}

	for (int i = 0; i < dim; i++){
		double min; config->getParameter("edu.uiuc.csl.system.var.min", i, &min);
		double max; config->getParameter("edu.uiuc.csl.system.var.max", i, &max);
		rrt.setBound(i, min, max);
	}

	for (int i = 0; i < param; i++){
		//rrt.setVariationBound(i, -0.005, 005);
		double min; config->getParameter("edu.uiuc.csl.system.param.min", i, &min);
		double max; config->getParameter("edu.uiuc.csl.system.param.max", i, &max);
		rrt.setVariationBound(i, min, max);
	}

	rrt.setConfig(config);
	rrt.setdt(dt);
	rrt.setSystem(circuit);

	cout << "config->get(\"edu.uiuc.crhc.core.mode\") = " << config->get("edu.uiuc.crhc.core.mode") << config->checkParameter("edu.uiuc.crhc.core.mode", "load") << endl;
	if (config->checkParameter("edu.uiuc.crhc.core.mode", "mc")){
		if (config->checkParameter("edu.uiuc.csl.core.sampling.digital", "1"))
			rrt.generateMonteCarloInputSequence();
		rrt.simulate(initialState);  //for just simulating the circuit
	}else if (config->checkParameter("edu.uiuc.crhc.core.mode", "rrt")){
		rrt.build(initialState);
	}else if (config->checkParameter("edu.uiuc.crhc.core.mode", "load")){
		rrt.load(config->get("edu.uiuc.crhc.core.options.inputFileName"));
	}else if (config->checkParameter("edu.uiuc.crhc.core.mode", "load+rrt")){
		rrt.load(config->get("edu.uiuc.crhc.core.options.inputFileName"));
		rrt.build();
	}
	else if (config->checkParameter("edu.uiuc.crhc.core.mode", "load+wca")){
		rrt.load(config->get("edu.uiuc.crhc.core.options.inputFileName"));
		rrt.worstCaseEyeDiagram();
	}
	else if (config->checkParameter("edu.uiuc.crhc.core.mode", "load+mc")){
		rrt.load(config->get("edu.uiuc.crhc.core.options.inputFileName"));
		rrt.simulate(iter, rrt.getRoot());
	}
	else if (config->checkParameter("edu.uiuc.crhc.core.mode", "mc+rrt")){
		rrt.simulate(initialState);
		rrt.build();
	}else{
		cout << "Uknown operation mode: " << config->get("edu.uiuc.crhc.core.mode") << endl;
	}


	cout << "Simulation finished " << endl;

	//rrt.save(outputFileName);

	//unit testing for the eye-diagram, this line can be safely commented or removed
	//rrt.getEyeDiagram()->test(); 
	
	if (config->checkParameter("edu.uiuc.crhc.core.options.plot", "1")){
		string plotPath; config->getParameter("edu.uiuc.crhc.core.options.plot.path", &plotPath);
		Plotter* plotter = new Plotter(plotPath);
		if (config->checkParameter("edu.uiuc.crhc.core.options.plot.type", "trace")){
			int v1; config->getParameter("edu.uiuc.crhc.core.options.plot.var[0]", &v1);
			string title = config->get("edu.uiuc.crhc.core.options.plot.title");
			plotter->plotTrace(rrt, v1, -1, dim , simulationTime, dt, title);
		}
		else if (config->checkParameter("edu.uiuc.crhc.core.options.plot.type", "rrt")){
			plotter->plotRRT("lines", rrt.getName().c_str(), "test", rrt, "x_1", "x_2", "t");
		}
		else if (config->checkParameter("edu.uiuc.crhc.core.options.plot.type", "eye")){

			stringstream str;
			double vmin = 0.7, vmax = 1.1, tmin = 0, tmax = 5e-10, voltage = 0, window = 0;
			double freq = 0;
			config->getParameter("edu.uiuc.csl.core.simulation.freq", &freq);
			config->getParameter("edu.uiuc.csl.core.simulation.window", &window);
			config->getParameter("edu.uiuc.crhc.core.options.eyediagram.var", &voltage);
			//config->getParameter("edu.uiuc.csl.system.var.min", voltage, &vmin);
			//config->getParameter("edu.uiuc.csl.system.var.max", voltage, &vmax);

			//tmax = window;
			double period = 1 / freq;


			str << "plot [ " << tmin << ":" << tmax << "][" << vmin << ":" << vmax << "] 0 with linespoints lt \"white\" pt 0.01";
			str << " title \"" << " " << "\"  \n";

			plotter->execute(str.str());

			//plotter->execute(rrt.getEyeDiagram()->toString());
			ofstream file;
			file.open("worst-eye-input-sequence");
			
			file << "Test Sequence for generating worst-case higher eyelid " << endl;
			vector<node*> fs = rrt.getEyeDiagram()->getFrontierSet(0);
			for (int i = 0; i < fs.size(); i++){
				vector<node*> path = rrt.getTest(fs[i]);
				for (int j = 0; j < path.size(); j++){
					file << i << " " << j << " " << path[j]->toString() << " ";
					for (int k = 0; k < path[j]->getInputVector().size(); k++){
						file << path[j]->getInput(k) << " ";
					}
					file << endl;
				}
				plotter->execute(rrt.drawTest(path, 1));
			}


			file << "Test Sequence for generating worst-case higher eyelid " << endl;
			vector<node*> fs2 = rrt.getEyeDiagram()->getFrontierSet(1);
			for (int i = 0; i < fs2.size(); i++){
				vector<node*> path = rrt.getTest(fs2[i]);
				for (int j = 0; j < path.size(); j++){
					file << i << " " << j << " " << path[j]->toString() << " ";

					for (int k = 0; k < path[j]->getInputVector().size(); k++){
						file << path[j]->getInput(k) << " ";
					}
					file << endl;
				}
				plotter->execute(rrt.drawTest(path, 2));
			}



			file.close();

		}
		else{
			cout << "Uknown plot command: [edu.uiuc.crhc.core.options.plot.type] " << config->get("edu.uiuc.crhc.core.options.plot.type") << endl;
		}
	}
	cout.flush();
	//rrt.getEyeDiagram()->sum();


	/*
	pll:
	vector<Monitor*> monitors;	//for this experiment, the monitors are empty.
	Plotter* plotter = new Plotter("C:\\opt\\gnuplot\\bin\\gnuplot.exe -persist");
	int mode = SIM_PLL;
	bool generatePlot = true;
	string inputFileName = "test2.rrt";
	string outputFileName = "pll_sim_ok_10000.rrt";
	//the new kernel has the monitors argument
	//kernel_RRT(monitors, mode, generatePlot, inputFileName, outputFileName, plotter);
	*/
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
	
	struct kdres *presults;
	//presults = kd_nearest_range(kd, q_sample->get(), 0.01); 
	//if(kd_res_size(presults)>0){	
	//	while( !kd_res_end( presults ) ) {
	//	// get the data and position of the current result item 
	//	node* res =  (node*) kd_res_item( presults, pos );
	//	cout << res->toString() << endl ;
	//	// compute the distance of the current result from the pt 
	//	dist = sqrt( dist_sq( pt, pos, 3 ) );

	//	// print out the retrieved data 
	//	printf( "node at (%.3f, %.3f, %.3f) is %.3f away and has data=%c\n", 
	//    pos[0], pos[1], pos[2], dist, *pch );

	//	// go to the next entry 
	//	kd_res_next( presults );
	//	}
	//   }else{
	 //       cout << "[error] 2" << endl ;
	//   }

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
	//kernel_RRT(monitors, mode, generatePlot, inputFileName, outputFileName, plotter);
	 

	//pt->parseFormula("");
	//pt->printParseTree(pt->getRoot());
	//setup an example tree execution to check the monitor
}


int main (int argc, const char * argv[]){
	srand((unsigned int)time(0));
	char full[_MAX_PATH];
	_fullpath(full, ".\\", _MAX_PATH);
	cout << "Current working directory is:" << full << endl;
	//string configFile = string(full) + "config\\half-wave-limiter.conf";
	//string configFile = string(full) + "config\\half-wave-limiter-81.conf";
	//string configFile = string(full) + "config\\josephson.conf";
	//string configFile = string(full) + "config\\inverter_mc_analysis.conf";
	//string configFile = string(full) + "config\\inverter_rrt_analysis.conf";
	string configFile = string(full) + "config\\tdo-ex1.conf";
	Configuration* config = new Configuration(configFile);

	//fft_experiments();
	//kdtree_experiment();
	//MonitorExperiment();
	kernel_RRT(config);
	//optimization_experiment();
	system("PAUSE");
	return 0;
}

