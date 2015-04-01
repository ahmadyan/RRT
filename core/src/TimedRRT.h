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

	vector<vector<node*> > nodeset; 
	vector<node*> frontier_set;
public:
	TimedRRT(Configuration* c, int _d, int _k, int, double _simTime, string);
	TimedRRT(Configuration* c, int, int, int, string);
	TimedRRT(Configuration* c, string fileName);
	node* findNearestNodeWithTimeIndex(node* q_sample, int, int);
	vector<double> generateSimulationParameters(node*, int);
	
	//generic methods
	void build(double* initialState);
	void build();
	void buildUniform();
	node* add_node(double*, int, int);

	//sample generation and search methods
	node* biasedSampling(double);

	void TimedRRT::simulate(double* initialState, vector< vector<double> > input_test);
	double getSimTime();
	void simulate(int iter, node* q_start, vector< vector<double> > test_input);
	void generateMonteCarloInputSequence();
	int worstCaseJitter(node* q_near);
	node* deltaSimulation(node* q_near);
	void worstCaseEyeDiagram();
	vector<double> loadPWLFile(string filename, int iterations, double);
	void compress_input();
	void construct_initial_frontier_set(int frontier_size, double frontier_set_min, double frontier_set_max);
	void update_frontier_set(node* q, int frontier_size, double frontier_set_min, double frontier_set_max);

	// --- model checking methods ---
	void addMonitor(Monitor* m);

	// --- test extraction methods ---
	vector<node*> computeGoalStates();
	vector<double> getTest(node*);
};
