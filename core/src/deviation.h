#pragma once

#include <iostream>
#include <vector>
#include "config.h"

using namespace std;

class Deviation{
	double dt;
	int size;
	vector< vector<double> > min; 
	vector< vector<double> > max; 
	public:
		Deviation(Configuration* , double);
		void reg();
		void push(double, double, int);
		double get(double time, int reg);
};