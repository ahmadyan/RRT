#pragma once

#include <vector>
#include <string>
#include "config.h"
#include "trace.h"
using namespace std;


enum SystemType {vanderpol, tunneldiode};
enum functionSign { positive, negative, posnegative, undefined };
double random(double, double);
int vanderpol_func (double t, const double y[], double f[], void *params);
int vanderpol_jac (double t, const double y[], double *dfdy, double dfdt[], void *params);
int tunneldiode_func (double t, const double y[], double f[], void *params);
int tunneldiode_jac  (double t, const double y[], double *dfdy, double dfdt[], void *params);

class Circuit{
    int dim ;
    double* min;
    double* max;
    double* init;
	double mu;
    SystemType type ;
public:
    
	Circuit(SystemType, Configuration*);
	~Circuit();
	
	double getMin(int);
	double getMax(int);
	double* getMinn();
	double* getMaxx();
	int getDimension();
	double* getInitialState();
	
	double* getVector(double* y);
	std::string generateVectorField();
    
	//transient simulation method, for a single state
    Trace* backward_simulate(Configuration*);
	Trace* simulate(Configuration*);
    double* forward_integration(double* , Configuration* config);
    double* backward_integration(double* , Configuration* config);
	//bool isReachable(Polytope* _source, Polytope* _sink);
	//functionSign reachability(Point* p1, Point* p2, int axis);
	//functionSign reachabilityUsingSampling(Point* p1, Point* p2);
	//functionSign fx(double x, double min, double max);
	//functionSign fy(double x, double min, double max);
};