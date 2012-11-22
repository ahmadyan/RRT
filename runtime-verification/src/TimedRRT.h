#pragma once
#include "Rrt.h"
#include "Monitor.h"
class TimedRRT : public RRT {
    double sim_time ;
	vector<Monitor*> monitors ;
public:
    TimedRRT(int _d, int _k, double _simTime, string);
    TimedRRT(int, int, string);
    TimedRRT(string fileName);
	void addMonitor(Monitor* m);
    void build(double* initialState, double);
	node* getNearestNode(node* q_sample);
	double getSimTime();
};