#pragma once
#include "Rrt.h"
#include "Monitor.h"
class TimedRRT : public RRT {
    double sim_time ;
	vector<Monitor*> monitors ;
public:
    TimedRRT(int _d, int _k, int , double _simTime, string);
    TimedRRT(int, int, int, string);
    TimedRRT(string fileName);
	vector<double> generateSimulationParameters(node*);
	void addMonitor(Monitor* m);
    void build(double* initialState);
	void simulate(double* initialState );
	double getSimTime();
};
