#include "vanderpol.h"
#include <fstream>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_ieee_utils.h>

Vanderpol::Vanderpol(Configuration* config){
	config->getParameter("edu.uiuc.csl.system.dimension", &d);
}

//Func will define the dynamics of the system. This is Van der Pol Oscillator
int vanderpol_func(double t, const double y[], double f[], void *params){
	vector<double> p = *(vector<double>*)params;
	double p0 = p.at(0);   
	double p1 = p.at(1);   
	double mu = 2;
	f[0] = y[1] + p0;
	f[1] = -(y[0] + p1) - mu*y[1] * (y[0] * y[0] - 1);
	return GSL_SUCCESS;
}

int vanderpol_jac(double t, const double y[], double *dfdy, double dfdt[], void *params){
	vector<double> p = *(vector<double>*)params;
	double p0 = p.at(0);   
	double p1 = p.at(1);   
	double mu = 2;
	gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 2, 2);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set(m, 0, 0, 0.0);
	gsl_matrix_set(m, 0, 1, 1.0);
	gsl_matrix_set(m, 1, 0, -2.0*mu*y[0] * y[1] - 1.0);
	gsl_matrix_set(m, 1, 1, -mu*(y[0] * y[0] - 1.0));
	dfdt[0] = 0.0;
	dfdt[1] = 0.0;
	return GSL_SUCCESS;
}

vector<double>  Vanderpol::simulate(double* ic, vector<double> param, vector<string> setting, double t0, double dt){
	vector<double> result;
	double t = t0;
	double tf = t0 + dt;
	result.push_back(tf);
	gsl_odeiv2_system sys = { vanderpol_func, vanderpol_jac, 2, &param };
	double y[2];
	for (int i = 0; i<2; i++) 
		y[i] = ic[i];

	double hstart = 1e-10;   //initial step size
	double epsabs = 1e-7;     //absolute error
	double epsrel = 1e-7;        //relative error

	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);
	//forward integration is usually easy:
	//— Function: int gsl_odeiv2_driver_apply (gsl_odeiv2_driver * d, double * t, const double t1, double y[])
	//  This function evolves the driver system d from t to t1. Initially vector y should contain the values of dependent variables at point t. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply is returned, and t and y contain the values from last successful step.
	//If maximum number of steps is reached, a value of GSL_EMAXITER is returned. If the step size drops below minimum value, the function returns with GSL_ENOPROG. If the user-supplied functions defined in the system sys returns GSL_EBADFUNC, the function returns immediately with the same return code. In this case the user must call gsl_odeiv2_driver_reset before calling this function again.

	int status = gsl_odeiv2_driver_apply(d, &t, tf, y);
	if (status != GSL_SUCCESS){
		printf("error, return value=%d\n", status);
	}
	gsl_odeiv2_driver_free(d);
	for (int i = 0; i < 2; i++){
		result.push_back(y[i]);
	}

	return result;

}

