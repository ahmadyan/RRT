#pragma once
#include <vector>
#include <iostream>
#include "hspiceInterface.h"
//#include <gsl/gsl_errno.h>
//#include <gsl/gsl_matrix.h>
//#include <gsl/gsl_odeiv2.h>

using namespace std;
class System{
    int d;  //dimension
    int (*function)(double t, const double y[], double f[], void *params);
    int (*jacobian) (double t, const double y[], double *dfdy, double dfdt[], void *params);
    CircuitType type;
public:
    System(CircuitType);
    ~System();
    System();
    
    void setSystem(int _d,int (*f)(double t, const double y[], double f[], void *params), 
                                int (*j) (double t, const double y[], double *dfdy, double dfdt[], void *params) );
    void dumpSystemInfo();
    double simulate();
    
    void check();
    void setDimension(int);
    int getDimension();
    double simulate(double* initialState, double* param, double);
	double integareODE(double* initialState, double* param, double);
	double runSPICE(double* initialState, double* param, double);
};
