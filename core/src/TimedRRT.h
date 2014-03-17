#pragma once
#include "Rrt.h"
#include "Monitor.h"
class TimedRRT : public RRT {
	double sim_time ;
	vector<Monitor*> monitors ;

	int* bits;
	double* jitter;
	double* transition;

	int GammaSimMode;
	double GammaJitter;
	double GammaTransition;
public:
	TimedRRT(Configuration* c, int _d, int _k, int, double _simTime, string);
	TimedRRT(Configuration* c, int, int, int, string);
	TimedRRT(Configuration* c, string fileName);
	vector<double> generateSimulationParameters(node*);
	void addMonitor(Monitor* m);
    void build(double* initialState);
	void simulate(double* initialState );
	double getSimTime();
	void build();
	void simulate(int iter, node* q_start);
	void generateMonteCarloInputSequence();
	int worstCaseJitter(node* q_near);
	void deltaSimulation(node* q_near);
	void worstCaseEyeDiagram();
};
