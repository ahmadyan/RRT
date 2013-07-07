#include "circuit.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <sstream>
#include <math.h>
#include <typeinfo>

#include "point.h"
#include "polytope.h"
#include "geometryUtil.h"
#include "rootFinding.h"

using namespace std;
using namespace geometry;


Circuit::Circuit(SystemType _type, Configuration* config){
    type=_type;
    dim  =   2;
    min = new double[dim];
    max = new double[dim];
    init= new double[dim];
    config->getParameter("edu.uiuc.crhc.core.system.x.min", &min[0]);
    config->getParameter("edu.uiuc.crhc.core.system.x.max", &max[0]);
    config->getParameter("edu.uiuc.crhc.core.system.y.min", &min[1]);
    config->getParameter("edu.uiuc.crhc.core.system.y.max", &max[1]);
    config->getParameter("edu.uiuc.crhc.core.system.x.init", &init[0]);
    config->getParameter("edu.uiuc.crhc.core.system.y.init", &init[1]);
    mu=2;
    if(type==vanderpol){
    }else if(type==tunneldiode){
    }
}
Circuit::~Circuit(){
    delete min;
    delete max;
    delete init;
}

double Circuit::getMin(int i){
    if(i<dim)
        return min[i];
    else
        return 0;
}

double Circuit::getMax(int i){
    if(i<dim)
        return max[i];
    else
        return 0;
}

int Circuit::getDimension(){
    return dim;
}

double* Circuit::getMinn(){
    return min;
}

double* Circuit::getMaxx(){
    return max ;
}

double* Circuit::getInitialState(){
    return init;
}

double* Circuit::getVector(double* y){
    double* f = new double[dim];
    for(int i=0;i<dim;i++) f[i]=0;
    if(type==vanderpol){
        vanderpol_func(0, y, f, &mu);
    }else if(type==tunneldiode){
        tunneldiode_func(0, y, f, &mu);
        
    }else{
        cout << "Circuit type is undefined" << endl ;
    }
    return f;
}

string Circuit::generateVectorField(){
    stringstream str;
    vector<double> p;
    p.push_back(0);
    p.push_back(0);
    double coef = 20;
    double dmax = 1 ;
    double* f = new double[dim];
    for(int i=0;i<dim;i++) f[i]=0;
    
    double* z = new double[dim];
    for(int i=0;i<dim;i++) z[i]=0;
    str << "set style line 1 lt 1 lw 3 pt 3 linecolor rgb \"blue\" " << endl ;
    str << "plot '-' using 1:2:3:4 ti \"phase portrait\" with vectors lt 2 head filled " << endl ;
    int pointScale=2;
    for(int i=1; i<pointScale*max[0]; i++){
        for(int j=1; j<pointScale*max[1]; j++){
            double x = (double)(i+min[0]); ///(double)pointScale ;
            double y = (double)(j+min[1]); ///(double)pointScale ;
            z[0] = x;
            z[1] = y;
            (type==vanderpol)? vanderpol_func(0, z, f, &p): tunneldiode_func(0, z, f, &p);
            double dx = f[0]/coef ;
            double dy = f[1]/coef ;
            
            if(dx>dmax) dx=dmax;
            if(dy>dmax) dy=dmax;
            if(dy<-dmax) dy=-dmax;
            if(dx<-dmax) dx=-dmax;
            
            str << x << " " << y << " " << dx << " " << dy << endl ;
        }
    }
    str << "e" << endl ;
    delete f;
    delete z;
    return str.str();
}

//------------------------------------------------------------
//
//    Transient Simulator for the System
//
//------------------------------------------------------------

double random(double a, double b){
    return (b-a)*(rand() / double(RAND_MAX)) + a;
}

//Func will define the dynamics of the system. This is Van der Pol Oscillator
int vanderpol_func (double t, const double y[], double f[], void *params){
    vector<double> p = *(vector<double>*)params;
    double p0= p.at(0);    //variation in diode
    double p1= p.at(1);    //variation in voltage source
    double mu=2;
    f[0] = y[1] + p0;
    f[1] = -(y[0]+p1) - mu*y[1]*(y[0]*y[0] - 1);
    return GSL_SUCCESS;
}

int vanderpol_jac (double t, const double y[], double *dfdy, double dfdt[], void *params){
    vector<double> p = *(vector<double>*)params;
    double p0= p.at(0);    //variation in diode
    double p1= p.at(1);    //variation in voltage source
    double mu=2;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    gsl_matrix * m = &dfdy_mat.matrix;
    gsl_matrix_set (m, 0, 0, 0.0);
    gsl_matrix_set (m, 0, 1, 1.0);
    gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    dfdt[0] = 0.0;
    dfdt[1] = 0.0;
    return GSL_SUCCESS;
}

int tunneldiode_func (double t, const double y[], double f[], void *params){
    vector<double> p = *(vector<double>*)params;
    double p0= p.at(0);    //variation in diode
    double p1= p.at(1);    //variation in voltage source
    double C=2;    //pF
    double L=5;    //uH
    double R=1.5;  //kOhm
    double u=1.2;  //V
    double vr = y[0];
    double il = y[1];
    double h = 17.76*vr -103.79*pow(vr,2) +229.62*pow(vr,3) -226.31*pow(vr,4) +83.72*pow(vr,5);
    f[0] = (1.0/C) * ( il - ( h + p0 ) ) ;
    f[1] = (1.0/L) * (-vr - R*il + u + p1 );
    return GSL_SUCCESS;
}

int tunneldiode_jac (double t, const double y[], double *dfdy, double dfdt[], void *params){
    cout << "implement me! [tunneldiode_jac]" << endl;
    double mu = 1 ;//*(double *)params;
    //gsl_matrix_view dfdy_mat = gsl_matrix_view_array (dfdy, 2, 2);
    //gsl_matrix * m = &dfdy_mat.matrix;
    //gsl_matrix_set (m, 0, 0, 0.0);
    //gsl_matrix_set (m, 0, 1, 1.0);
    //gsl_matrix_set (m, 1, 0, -2.0*mu*y[0]*y[1] - 1.0);
    //gsl_matrix_set (m, 1, 1, -mu*(y[0]*y[0] - 1.0));
    //dfdt[0] = 0.0;
    //dfdt[1] = 0.0;
    //    return GSL_SUCCESS;
    return -2;
}

Trace* Circuit::simulate(Configuration* config){
    Trace* trace = new Trace();
    
    vector<double> p;
    p.push_back(0);
    p.push_back(0);
    
    
    gsl_odeiv2_system sys;
    if(type==vanderpol) {
        sys = {vanderpol_func, vanderpol_jac, 2, &p};
    }else if(type==tunneldiode){
        sys = {tunneldiode_func, tunneldiode_jac, 2, &p};
    }else{
        cout << "undefined circuit" << endl ;
    }
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-9, 1e-9, 0.0);
    
    int i;
    double t = 0.0, t1 = 100.0;
    //config->getParameter("edu.uiuc.crhc.core.simulation.totalTime", &t1);
    double* state = getInitialState();
    double y[2] = { state[0], state[1] };
    
    trace->AddSample(make_pair(y[0], y[1]));
    
    for (int i = 1; i <= 1000; i++){
        double ti = i * t1 / 1000.0;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        trace->AddSample(make_pair(y[0], y[1]));
    }
    gsl_odeiv2_driver_free (d);
    return trace;
}

double* Circuit::forward_integration(double* state, Configuration* config){
    vector<double> p;
    double var0min=-1;
    double var0max=1;
    double var1min=-1;
    double var1max=1;
    config->getParameter("edu.uiuc.crhc.core.system.variation.0.min", &var0min);
    config->getParameter("edu.uiuc.crhc.core.system.variation.0.max", &var0max);
    config->getParameter("edu.uiuc.crhc.core.system.variation.1.min", &var1min);
    config->getParameter("edu.uiuc.crhc.core.system.variation.1.max", &var1max);
    
    double p0 = random(var0min, var0max);
    double p1 = random(var1min, var1max);
    p.push_back(p0);
    p.push_back(p1);
    
    gsl_odeiv2_system sys;
    if(type==vanderpol) {
        sys = {vanderpol_func, vanderpol_jac, 2, &p};
    }else if(type==tunneldiode){
        sys = {tunneldiode_func, NULL, 2, &p};
    }else{
        cout << "undefined circuit" << endl ;
    }
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-9, 1e-9, 0.0);
    double t = 0.0, tf = 1;
    config->getParameter("edu.uiuc.crhc.core.simulation.step", &tf);
    double y[2] = { state[0], state[1] };
    int status = gsl_odeiv2_driver_apply (d, &t, tf, y);
    //int status =  gsl_odeiv2_driver_apply_fixed_step (d,&t, tf, 1, y);
    
    if (status != GSL_SUCCESS){
        printf ("error, return value=%d\n", status);
        cout << y[0] << " " << y[1] << "  --> " << state[0] << " " << state[1] << endl ;
    }
    gsl_odeiv2_driver_free (d);
    
    double* result = new double[2];
    for(int i=0;i<2;i++){
        result[i] = y[i];
    }
    //cout << "resulting points are " << result[0] << " " << result[1] << endl ;
    return result ;
}

double* Circuit::backward_integration(double* state, Configuration* config){
    vector<double> p;
    double var0min=-1;
    double var0max=1;
    double var1min=-1;
    double var1max=1;
    config->getParameter("edu.uiuc.crhc.core.system.variation.0.min", &var0min);
    config->getParameter("edu.uiuc.crhc.core.system.variation.0.max", &var0max);
    config->getParameter("edu.uiuc.crhc.core.system.variation.1.min", &var1min);
    config->getParameter("edu.uiuc.crhc.core.system.variation.1.max", &var1max);
    
    double p0 = random(var0min, var0max);
    double p1 = random(var1min, var1max);
    
    p.push_back(p0);
    p.push_back(p1);
    gsl_odeiv2_system sys;
    if(type==vanderpol) {
        sys = {vanderpol_func, vanderpol_jac, 2, &p};
    }else if(type==tunneldiode){
        sys = {tunneldiode_func, NULL, 2, &p};
    }else{
        cout << "undefined circuit" << endl ;
    }
    
    double t = 0.0, tf = 1;
    config->getParameter("edu.uiuc.crhc.core.simulation.backwardstep", &tf);
    //cout << "simulating circuit from " << y[0] << " " << y[1] << " time t0= " << t << " for " << tf << endl ;
    //int status = gsl_odeiv2_driver_apply (d, &t, -tf, y);
    
    
    bool isNan=true ;
    int status=0;
    double* result = new double[2];
    while(isNan){
        gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-9, 1e-9, 0.0);
        double y[2] = { state[0], state[1] };
        status =  gsl_odeiv2_driver_apply_fixed_step (d,&t, -tf, 1, y);
        isNan =  isnan(y[0]) || isnan(y[1]);
        if(isNan) cout << "reducing integration time" << endl ;
        tf = tf/10.0;
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            isNan = true;
            cout << "reducing integration time" << endl ;
            cout << y[0] << " " << y[1] << "  --> " << state[0] << " " << state[1] << endl ;
        }
        for(int i=0;i<2;i++){
            result[i] = y[i];
        }
        gsl_odeiv2_driver_free (d);
    }
    
    
    
    
    //cout << "resulting points are " << result[0] << " " << result[1] << endl ;
    return result ;
}

Trace* Circuit::backward_simulate(Configuration* config){
    Trace* trace = new Trace();
    gsl_odeiv2_system sys;
    if(type==vanderpol) {
        sys = {vanderpol_func, vanderpol_jac, 2, &mu};
    }else if(type==tunneldiode){
        sys = {tunneldiode_func, tunneldiode_jac, 2, &mu};
    }else{
        cout << "undefined circuit" << endl ;
    }
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double t = 0.0, t1 = 0;
    config->getParameter("edu.uiuc.crhc.core.simulation.backwardstep", &t1);
    double* state = getInitialState();
    double y[2] = { state[0], state[1] };
    for (int i = 1; i <= 1000; i++){
        double ti = i * t1 / 1000.0;
        /*
         — Function: int gsl_odeiv2_driver_apply (gsl_odeiv2_driver * d, double * t, const double t1, double y[])
         This function evolves the driver system d from t to t1. Initially vector y should contain the values of dependent variables at point t. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply is returned, and t and y contain the values from last successful step.
         
         If maximum number of steps is reached, a value of GSL_EMAXITER is returned. If the step size drops below minimum value, the function returns with GSL_ENOPROG. If the user-supplied functions defined in the system sys returns GSL_EBADFUNC, the function returns immediately with the same return code. In this case the user must call gsl_odeiv2_driver_reset before calling this function again.
         
         — Function: int gsl_odeiv2_driver_apply_fixed_step (gsl_odeiv2_driver * d, double * t, const double h, const unsigned long int n, double y[])
         This function evolves the driver system d from t with n steps of size h. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply_fixed_step is returned, and t and y contain the values from last successful step.
         */
        
        //int status = gsl_odeiv2_driver_apply (d, &t, -ti, y);
        int status =  gsl_odeiv2_driver_apply_fixed_step (d,&t, -ti, 1, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        trace->AddSample(make_pair(y[0], y[1]));
    }
    gsl_odeiv2_driver_free (d);
    return trace;
}