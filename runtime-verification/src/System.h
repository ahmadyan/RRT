#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <vector>
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_odeiv2.h>

using namespace std;
enum simulatorType {SPICE, GSL};
class System{
	
public:
	~System();
	System();

	simulatorType simulator;
	int d;	//dimension, currently not in use since the ode integration id depreted
	int(*function)(double t, const double y[], double f[], void *params);
	int(*jacobian) (double t, const double y[], double *dfdy, double dfdt[], void *params);

	static vector<double> parse(string str);
	static bool is_only_ascii_whitespace(const std::string& str);
	static double unit(char u);

	///	Simulate is the most important function in the system class. It should be called from the child class (like PLL, TDO or Inv)
	/// In the input, we read the initial condition as state, the list of inputs as parameters and other settings such as filename, etc.
	/// The return is asymetrical. The result[0] should contains the simulation time (usually dt, unless specified otherwise). Then result[i], i>0
	/// should contain the states of the simulation
	virtual vector<double> simulate(double* state, vector<double> parameters, vector<string> settings, double t0,  double dt) = 0;
	
    void setSystem(int _d,int (*f)(double t, const double y[], double f[], void *params), 
                          int (*j) (double t, const double y[], double *dfdy, double dfdt[], void *params) );
    
    void check();
    void setDimension(int);
    int getDimension();
    //double simulate(double* initialState, double* param, double);
	double integareODE(double* initialState, double* param, double);
	//double runSPICE(double* initialState, double* param, double);
};
