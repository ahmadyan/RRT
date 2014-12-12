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
#include "analogprop.h"

#include "System.h"
#include "spice.h"
#include "vanderpol.h"
#include "pll.h"
#include "tdo.h"
#include "inverter.h"
#include "josephson.h"

using namespace std;

System* systemSelector(Configuration* config){
	System* circuit;
	cout << config->get("edu.uiuc.csl.system.name") << endl;
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
	string outputFileName;  config->getParameter("edu.uiuc.csl.core.options.outputFileName", &outputFileName);
	string inputFileName;  config->getParameter("edu.uiuc.csl.core.options.inputFileName", &inputFileName);

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

	// if the mode of operation is Monte Carlo or Test compression, compute the number of iterations by dividing the simulation time over dt. 
	// In zip mode, the simulation time denote to the length of the test sequence that we wish to minimize
	int iter = 0;
	if (config->checkParameter("edu.uiuc.csl.core.mode", "mc") || config->checkParameter("edu.uiuc.csl.core.mode", "zip")){
		iter = simulationTime / dt;
	}else{
		//otherwise (if we are loading, using rrt, etc) directly get the number of iteration from the iteration variable.
		config->getParameter("edu.uiuc.csl.core.sampling.iteration", &iter);
	}

	cout << "----------------------------------------------" << endl;
	cout << "RRT Options: " << endl;
	cout << "Iter=" << iter << endl;
	cout << "sim-time=" << simulationTime << endl;
	cout << "dt=" << dt << endl;
	cout << "----------------------------------------------" << endl;

	//Construct the timed rrt data structure
	TimedRRT rrt = TimedRRT(config, dim, iter, param, simulationTime, name);

	//initialize the parameter space and it's boundaries
	for (int i = 0; i < param; i++){		
		double min; config->getParameter("edu.uiuc.csl.system.param.min", i, &min);
		double max; config->getParameter("edu.uiuc.csl.system.param.max", i, &max);
		rrt.setVariationBound(i, min, max);
	}

	for (int i = 0; i < dim; i++){
		double min; config->getParameter("edu.uiuc.csl.system.var.min", i, &min);
		double max; config->getParameter("edu.uiuc.csl.system.var.max", i, &max);
		rrt.setBound(i, min, max);
	}


	//runtime monitoring experiments
	//the model checking engine will be loaded with the edu.uiuc.csl.validation variable
	if (config->checkParameter("edu.uiuc.csl.validation.enable", "1")){
		vector<Monitor*> monitors;
		for (int i = 1; i <= 4; i++){//property loop, todo: implement this from a text parser. not hardcoded
			monitors.push_back(new Monitor(new Property(i)));
		}
		//AnalogProperty* ap = (AnalogProperty*)(monitors[0]->property->argument);
		rrt.setMonitor(monitors);

		//for (int i = 0; i<(int)(monitors.size()); i++){
		//	rrt.addMonitor(monitors[i]);
		//}
	}
	
	rrt.setConfig(config);
	rrt.setdt(dt);
	rrt.setSystem(circuit);

	cout << "config->get(\"edu.uiuc.csl.core.mode\") = " << config->get("edu.uiuc.csl.core.mode") << config->checkParameter("edu.uiuc.csl.core.mode", "load") << endl;
	if (config->checkParameter("edu.uiuc.csl.core.mode", "mc")){
		if (config->checkParameter("edu.uiuc.csl.core.sampling.digital", "1"))
			rrt.generateMonteCarloInputSequence();
		rrt.simulate(initialState, vector<vector<double>>());  //for just simulating the circuit
	}else if (config->checkParameter("edu.uiuc.csl.core.mode", "rrt")){
		rrt.build(initialState);
	}else if (config->checkParameter("edu.uiuc.csl.core.mode", "load")){
		rrt.load(config->get("edu.uiuc.csl.core.options.inputFileName"));
	}else if (config->checkParameter("edu.uiuc.csl.core.mode", "load+rrt")){
		rrt.load(config->get("edu.uiuc.csl.core.options.inputFileName"));
		rrt.build();
	}
	else if (config->checkParameter("edu.uiuc.csl.core.mode", "load+wca")){
		rrt.load(config->get("edu.uiuc.csl.core.options.inputFileName"));
		rrt.worstCaseEyeDiagram();
	}
	else if (config->checkParameter("edu.uiuc.csl.core.mode", "load+mc")){
		rrt.load(config->get("edu.uiuc.csl.core.options.inputFileName"));
		rrt.simulate(iter, rrt.getRoot(), vector<vector<double>>());
	}
	else if (config->checkParameter("edu.uiuc.csl.core.mode", "mc+rrt")){
		rrt.simulate(initialState, vector<vector<double>>());
		rrt.build();
	}
	else if (config->checkParameter("edu.uiuc.csl.core.mode", "zip")){
		cout << "Test compression engine..." << endl;
		//load initial test sequence
		vector< vector<double> > input_test;
		for (int i = 0; i < param; i++){
			if (config->checkParameter("edu.uiuc.csl.system.param.type", i, "file")){
				string filename; 
				config->getParameter("edu.uiuc.csl.system.param.file", i, &filename);
				vector<double> pwl_timeseries = rrt.loadPWLFile(filename, iter, dt);
				input_test.push_back(pwl_timeseries);
			}
		}
		
		rrt.simulate(initialState, input_test);
		rrt.compress_input();
	}
	else if (config->checkParameter("edu.uiuc.csl.core.mode", "load+zip")){
		rrt.load(config->get("edu.uiuc.csl.core.options.inputFileName"));
		rrt.compress_input();
	}else{
		cout << "Uknown operation mode: " << config->get("edu.uiuc.csl.core.mode") << endl;
	}


	cout << "Simulation finished " << endl;
	rrt.save(outputFileName);

	//unit testing for the eye-diagram, this line can be safely commented or removed
	//rrt.getEyeDiagram()->test(); 
	
	if (config->checkParameter("edu.uiuc.csl.core.options.plot", "1")){
		string plotPath; config->getParameter("edu.uiuc.csl.core.options.plot.path", &plotPath);
		Plotter* plotter = new Plotter(plotPath);
		if (config->checkParameter("edu.uiuc.csl.core.options.plot.type", "trace")){
			int v1; config->getParameter("edu.uiuc.csl.core.options.plot.var[0]", &v1);
			string title = config->get("edu.uiuc.csl.core.options.plot.title");
			plotter->plotTrace(rrt, v1, -1, dim , simulationTime, dt, title);
		}
		else if (config->checkParameter("edu.uiuc.csl.core.options.plot.type", "rrt")){
			plotter->plotRRT("lines", rrt.getName().c_str(), "test", rrt, "x_1", "x_2", "t");
		}
		else if (config->checkParameter("edu.uiuc.csl.core.options.plot.type", "eye")){

			stringstream str;
			double vmin = 0.7, vmax = 1.1, tmin = 0, tmax = 5e-10, voltage = 0, window = 0;
			double freq = 0;
			config->getParameter("edu.uiuc.csl.core.simulation.freq", &freq);
			config->getParameter("edu.uiuc.csl.core.simulation.window", &window);
			config->getParameter("edu.uiuc.csl.core.options.eyediagram.var", &voltage);
			//config->getParameter("edu.uiuc.csl.system.var.min", voltage, &vmin);
			//config->getParameter("edu.uiuc.csl.system.var.max", voltage, &vmax);

			//tmax = window;
			double period = 1 / freq;


			str << "plot [ " << tmin << ":" << tmax << "][" << vmin << ":" << vmax << "] 0 with linespoints lt \"white\" pt 0.01";
			str << " title \"" << " " << "\"  \n";

			//plotter->execute(str.str());

			//plotter->execute(rrt.getEyeDiagram()->toString());
			//plotter->execute(rrt.getEyeDiagram()->drawContour());

			/*
			// We used this piece of code for generating worst-case input perturbation sequence that resulted in the worst-case eye diagram
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
			*/

			

		}
		else{
			cout << "Uknown plot command: [edu.uiuc.csl.core.options.plot.type] " << config->get("edu.uiuc.csl.core.options.plot.type") << endl;
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


void fft_experiments(){
	//	Computing the joint time-frequency space instead of only time-augmented RRT
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

int main (int argc, const char * argv[]){
	srand((unsigned int)time(0));
	char full[_MAX_PATH];
	_fullpath(full, ".\\", _MAX_PATH);
	cout << "Current working directory is:" << full << endl;
	string configFile = string(full) + "config//Opamp_single_ended_unit_gain.conf";
	Configuration* config = new Configuration(configFile);
	kernel_RRT(config);
	system("PAUSE");
	return 0;
}
