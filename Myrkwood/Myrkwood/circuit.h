#pragma once

#include <vector>
#include <string>
//#include "polytope.h"
using namespace std;


enum SystemType {vanderpol, tunneldiode};
enum functionSign { positive, negative, posnegative, undefined };

int func (double t, const double y[], double f[], void *params);
int jac (double t, const double y[], double *dfdy, double dfdt[], void *params);

class Circuit{
	 int dim ;
     double* min ;
     double* max ;
	double mu;

public:

	Circuit(SystemType);
	~Circuit();
	
	double getMin(int);
	double getMax(int);
	double* getMinn();
	double* getMaxx();
	int getDimension();
	double* getInitialState();
	double random(double, double);
	double* getVector(double* y);
	std::string generateVectorField();
        
	//transient simulation method, for a single state
	std::vector<std::pair<double, double> > simulate();

	//bool isReachable(Polytope* _source, Polytope* _sink);
	//functionSign reachability(Point* p1, Point* p2, int axis);
	//functionSign reachabilityUsingSampling(Point* p1, Point* p2);
	//functionSign fx(double x, double min, double max);
	//functionSign fy(double x, double min, double max);
};