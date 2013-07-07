#include "circuit.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_ieee_utils.h>
#include <sstream>
#include <math.h>
#include <typeinfo>

//Radau for DAEs
#include "StiffIntegratorT.h"

#include "point.h"
#include "polytope.h"
#include "geometryUtil.h"
#include "rootFinding.h"

using namespace std;
using namespace geometry;


/* Maximum number of ODE equations */
#define MAXEQ 15

/* Maximum number of ODE solvers */
#define MAXNS 20

/* Track number of function and jacobian evaluations in tests
 with global variables
 */
int nfe;
int nje;

vector<double> variation ;
Circuit::Circuit(SystemType _type, Configuration* config){
    switch (_type) {
        case vanderpol:
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

            break;
        case tunneldiode:
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

            break;
        case ringmodulator:
            type=_type;
            dim  =   15;
            min = new double[dim];
            max = new double[dim];
            init= new double[dim];
            for(int i=0;i<dim;i++){
                min[i]=-1;
                max[i]=1;
                init[i]=0;
            }
            break;
        case dae_amplifier:
            type=_type;
            dim=8;
            min = new double[dim];
            max = new double[dim];
            init= new double[dim];
            
            for(int i=0;i<dim;i++){
                min[i]=-5;
                max[i]=5;
                init[i]=0;
            }
            init[0] = 0;
            init[1] = 6;
            init[2] = 3;
            init[3] = 3;
            init[4] = 6;
            init[5] = 3;
            init[6] = 3;
            init[7] = 0;

            break;
        default:
            break;
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


/* Ring Modulator, stiff ODE of dimension 15.
 
 Reference: Walter M. Lioen, Jacques J.B. de Swart, Test Set for
 Initial Value Problem Solvers, Release 2.1 September 1999,
 http://ftp.cwi.nl/IVPtestset/software.htm
 */

#define NRINGMOD 15
#define MAXD 15

void Function(double x, double *y, double *f)
{
	double ue, ub, uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
	double w, uet, fac1, fac2;
    
	ue = 0.1;
	ub = 6.0;
	uf = 0.026;
	alpha = 0.99;
	beta = 1.0e-6;
	r0 = 1000.0;
	r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = 9000.0;
	w = 2.0*3.141592654*100.0;
	uet = ue*sin(w*x) + variation[0];
	fac1 = beta*(exp((y[3] - y[2])/uf) - 1.0);
	fac2 = beta*(exp((y[6] - y[5])/uf) - 1.0);
    
	f[0] = y[0]/r9;
	f[1] = (y[1] - ub)/r8 + alpha*fac1;
	f[2] = y[2]/r7 - fac1;
	f[3] = y[3]/r5 + (y[3] - ub)/r6 + (1.0 - alpha)*fac1;
	f[4] = (y[4] - ub)/r4 + alpha*fac2;
	f[5] = y[5]/r3 - fac2;
	f[6] = y[6]/r1 + (y[6] - ub)/r2 + (1.0 - alpha)*fac2;
	f[7] = (y[7] - uet)/r0;
    
	return;
    
} // Function

void Jacobian(double x, double *y, double **J)
{
    
	double uf, alpha, beta, r0, r1, r2, r3, r4, r5, r6, r7, r8, r9;
	double fac14, fac27;
    
	uf = 0.026;
	alpha = 0.99;
	beta = 1.0e-6;
	r0 = 1000.0;
	r1 = r2 = r3 = r4 = r5 = r6 = r7 = r8 = r9 = 9000.0;
	fac14 = beta*exp((y[3] - y[2])/uf)/uf;
	fac27 = beta*exp((y[6] - y[5])/uf)/uf;
    
	for (int j = 0; j < 8; j++) {
		J[0][j] = 0.0;
		J[1][j] = 0.0;
		J[3][j] = 0.0;
	}
	J[2][0] = 1.0/r9;
	J[2][1] = 1.0/r8;
	J[1][2] = -alpha*fac14;
	J[0][3] = alpha*fac14;
	J[2][2] = 1.0/r7 + fac14;
	J[1][3] = -fac14;
	J[2][3] = 1.0/r5 + 1.0/r6 + (1.0 - alpha)*fac14;
	J[3][2] = -(1.0 - alpha)*fac14;
	J[2][4] = 1.0/r4;
	J[1][5] = -alpha*fac27;
	J[0][6] = alpha*fac27;
	J[2][5] = 1.0/r3 + fac27;
	J[1][6] = -fac27;
	J[2][6] = 1.0/r1 + 1.0/r2 + (1.0 - alpha)*fac27;
	J[3][5] = -(1.0 - alpha)*fac27;
	J[2][7] = 1.0/r0;
    
	return;
    
} // Jacobian

void Mass(double **M)
{
 	const double c1 = 1.0e-6;
	const double c2 = 2.0e-6;
	const double c3 = 3.0e-6;
	const double c4 = 4.0e-6;
	const double c5 = 5.0e-6;
    
	for (int i = 0; i < 8; i++) {
		M[0][i] = 0.0;
		M[2][i] = 0.0;
	}
    
	M[1][0] = -c5;
	M[0][1] = c5;
	M[2][0] = c5;
	M[1][1] = -c5;
	M[1][2] = -c4;
	M[1][3] = -c3;
	M[0][4] = c3;
	M[2][3] = c3;
	M[1][4] = -c3;
	M[1][5] = -c2;
	M[1][6] = -c1;
	M[0][7] = c1;
	M[2][6] = c1;
	M[1][7] = -c1;
    
	return;
    
} // Mass


int rhs_ringmod (double t, const double y[], double f[], void *params){
    vector<double> p = *(vector<double>*)params;
    const double c = 1.6e-8;
    const double cs = 2e-12;
    const double cp = 1e-8;
    const double r = 25e3;
    const double rp = 50e0;
    const double lh = 4.45e0;
    const double ls1 = 2e-3;
    const double ls2 = 5e-4;
    const double ls3 = 5e-4;
    const double rg1 = 36.3;
    const double rg2 = 17.3;
    const double rg3 = 17.3;
    const double ri = 5e1;
    const double rc = 6e2;
    const double gamma = 40.67286402e-9;
    const double delta = 17.7493332;
    const double pi = 3.141592653589793238462643383;
    
    const double uin1 = 0.5 * sin (2e3 * pi * t) + p.at(0);
    const double uin2 = 2 * sin (2e4 * pi * t) + p.at(1);
    //const double uin2 = sin (2e4 * pi * t) + p.at(1);
    const double ud1 = +y[2] - y[4] - y[6] - uin2 ;// + p.at(2);
    const double ud2 = -y[3] + y[5] - y[6] - uin2 ;//+ p.at(3);
    const double ud3 = +y[3] + y[4] + y[6] + uin2 ;//+ p.at(4);
    const double ud4 = -y[2] - y[5] + y[6] + uin2 ;//+ p.at(5);
    
    const double qud1 = gamma * (exp (delta * ud1) - 1.0);
    const double qud2 = gamma * (exp (delta * ud2) - 1.0);
    const double qud3 = gamma * (exp (delta * ud3) - 1.0);
    const double qud4 = gamma * (exp (delta * ud4) - 1.0);
    
    extern int nfe;
    nfe += 1;
    
    f[0] = (y[7] - 0.5 * y[9] + 0.5 * y[10] + y[13] - y[0] / r) / c;
    f[1] = (y[8] - 0.5 * y[11] + 0.5 * y[12] + y[14] - y[1] / r) / c;
    f[2] = (y[9] - qud1 + qud4) / cs;
    f[3] = (-y[10] + qud2 - qud3) / cs;
    f[4] = (y[11] + qud1 - qud3) / cs;
    f[5] = (-y[12] - qud2 + qud4) / cs;
    f[6] = (-y[6] / rp + qud1 + qud2 - qud3 - qud4) / cp;
    f[7] = -y[0] / lh;
    f[8] = -y[1] / lh;
    f[9] = (0.5 * y[0] - y[2] - rg2 * y[9]) / ls2;
    f[10] = (-0.5 * y[0] + y[3] - rg3 * y[10]) / ls3;
    f[11] = (0.5 * y[1] - y[4] - rg2 * y[11]) / ls2;
    f[12] = (-0.5 * y[1] + y[5] - rg3 * y[12]) / ls3;
    f[13] = (-y[0] + uin1 - (ri + rg1) * y[13]) / ls1;
    f[14] = (-y[1] - (rc + rg1) * y[14]) / ls1;
    
    return GSL_SUCCESS;
}

int jac_ringmod (double t, const double y[], double *dfdy, double dfdt[], void *params){
    vector<double> p = *(vector<double>*)params;
    const double c = 1.6e-8;
    const double cs = 2e-12;
    const double cp = 1e-8;
    const double r = 25e3;
    const double rp = 50e0;
    const double lh = 4.45e0;
    const double ls1 = 2e-3;
    const double ls2 = 5e-4;
    const double ls3 = 5e-4;
    const double rg1 = 36.3;
    const double rg2 = 17.3;
    const double rg3 = 17.3;
    const double ri = 5e1;
    const double rc = 6e2;
    const double gamma = 40.67286402e-9;
    const double delta = 17.7493332;
    const double pi = 3.141592653589793238462643383;
    
    const double uin2 = 2 * sin (2e4 * pi * t) + p.at(1);
    //const double uin2 =  sin (2e4 * pi * t) + p.at(1);
    const double ud1 = +y[2] - y[4] - y[6] - uin2 ;//+ p.at(2);
    const double ud2 = -y[3] + y[5] - y[6] - uin2 ;//+ p.at(3);
    const double ud3 = +y[3] + y[4] + y[6] + uin2 ;//+ p.at(4);
    const double ud4 = -y[2] - y[5] + y[6] + uin2 ;//+ p.at(5);
    const double qpud1 = gamma * delta * exp (delta * ud1);
    const double qpud2 = gamma * delta * exp (delta * ud2);
    const double qpud3 = gamma * delta * exp (delta * ud3);
    const double qpud4 = gamma * delta * exp (delta * ud4);
    
    extern int nje;
    size_t i;
    
    nje += 1;
    
    for (i = 0; i < NRINGMOD * NRINGMOD; i++)
    {
        dfdy[i] = 0.0;
    }
    
    
    
    dfdy[0 * NRINGMOD + 0] = -1 / (c * r);
    dfdy[0 * NRINGMOD + 7] = 1 / c;
    dfdy[0 * NRINGMOD + 9] = -0.5 / c;
    dfdy[0 * NRINGMOD + 10] = -dfdy[0 * NRINGMOD + 9];
    dfdy[0 * NRINGMOD + 13] = dfdy[0 * NRINGMOD + 7];
    dfdy[1 * NRINGMOD + 1] = dfdy[0 * NRINGMOD + 0];
    dfdy[1 * NRINGMOD + 8] = dfdy[0 * NRINGMOD + 7];
    dfdy[1 * NRINGMOD + 11] = dfdy[0 * NRINGMOD + 9];
    dfdy[1 * NRINGMOD + 12] = dfdy[0 * NRINGMOD + 10];
    dfdy[1 * NRINGMOD + 14] = dfdy[0 * NRINGMOD + 13];
    dfdy[2 * NRINGMOD + 2] = (-qpud1 - qpud4) / cs;
    dfdy[2 * NRINGMOD + 4] = qpud1 / cs;
    dfdy[2 * NRINGMOD + 5] = -qpud4 / cs;
    dfdy[2 * NRINGMOD + 6] = (qpud1 + qpud4) / cs;
    dfdy[2 * NRINGMOD + 9] = 1 / cs;
    dfdy[3 * NRINGMOD + 3] = (-qpud2 - qpud3) / cs;
    dfdy[3 * NRINGMOD + 4] = -qpud3 / cs;
    dfdy[3 * NRINGMOD + 5] = qpud2 / cs;
    dfdy[3 * NRINGMOD + 6] = (-qpud2 - qpud3) / cs;
    dfdy[3 * NRINGMOD + 10] = -1 / cs;
    dfdy[4 * NRINGMOD + 2] = qpud1 / cs;
    dfdy[4 * NRINGMOD + 3] = -qpud3 / cs;
    dfdy[4 * NRINGMOD + 4] = (-qpud1 - qpud3) / cs;
    dfdy[4 * NRINGMOD + 6] = (-qpud1 - qpud3) / cs;
    dfdy[4 * NRINGMOD + 11] = 1 / cs;
    dfdy[5 * NRINGMOD + 2] = -qpud4 / cs;
    dfdy[5 * NRINGMOD + 3] = qpud2 / cs;
    dfdy[5 * NRINGMOD + 5] = (-qpud2 - qpud4) / cs;
    dfdy[5 * NRINGMOD + 6] = (qpud2 + qpud4) / cs;
    dfdy[5 * NRINGMOD + 12] = -1 / cs;
    dfdy[6 * NRINGMOD + 2] = (qpud1 + qpud4) / cp;
    dfdy[6 * NRINGMOD + 3] = (-qpud2 - qpud3) / cp;
    dfdy[6 * NRINGMOD + 4] = (-qpud1 - qpud3) / cp;
    dfdy[6 * NRINGMOD + 5] = (qpud2 + qpud4) / cp;
    dfdy[6 * NRINGMOD + 6] = (-qpud1 - qpud2 - qpud3 - qpud4 - 1 / rp) / cp;
    dfdy[7 * NRINGMOD + 0] = -1 / lh;
    dfdy[8 * NRINGMOD + 1] = dfdy[7 * NRINGMOD + 0];
    dfdy[9 * NRINGMOD + 0] = 0.5 / ls2;
    dfdy[9 * NRINGMOD + 2] = -1 / ls2;
    dfdy[9 * NRINGMOD + 9] = -rg2 / ls2;
    dfdy[10 * NRINGMOD + 0] = -0.5 / ls3;
    dfdy[10 * NRINGMOD + 3] = 1 / ls3;
    dfdy[10 * NRINGMOD + 10] = -rg3 / ls3;
    dfdy[11 * NRINGMOD + 1] = dfdy[9 * NRINGMOD + 0];
    dfdy[11 * NRINGMOD + 4] = dfdy[9 * NRINGMOD + 2];
    dfdy[11 * NRINGMOD + 11] = dfdy[9 * NRINGMOD + 9];
    dfdy[12 * NRINGMOD + 1] = dfdy[10 * NRINGMOD + 0];
    dfdy[12 * NRINGMOD + 5] = dfdy[10 * NRINGMOD + 3];
    dfdy[12 * NRINGMOD + 12] = dfdy[10 * NRINGMOD + 10];
    dfdy[13 * NRINGMOD + 0] = -1 / ls1;
    dfdy[13 * NRINGMOD + 13] = -(ri + rg1) / ls1;
    dfdy[14 * NRINGMOD + 1] = dfdy[13 * NRINGMOD + 0];
    dfdy[14 * NRINGMOD + 14] = -(rc + rg1) / ls1;
    
    for (i = 0; i < NRINGMOD; i++)
    {
        dfdt[i] = 0.0;
    }
    
    return GSL_SUCCESS;
}

// For nonstiff system, we use rk8pd stepper with a fixed step size.
// This function works for time-invariant systems (t=0)
double* Circuit::integration_nonstiff(gsl_odeiv2_system sys, double* state, Configuration* config, bool direction, double t0){
    double t = t0;
    double tf = 1e-14;
    if(direction){
        config->getParameter("edu.uiuc.crhc.core.simulation.step", &tf);
    }else{
        config->getParameter("edu.uiuc.crhc.core.simulation.backwardstep", &tf);
    }
    tf += t0;
   
    double y[MAXEQ];
    for(int i=0;i<MAXEQ;i++) y[i]=state[i];
    
    if(direction){
        for(int i=0;i<MAXEQ;i++) y[i]=state[i];
        double hstart=1e-10;   //initial step size
        double epsabs=1e-7;     //absolute error
        double epsrel=1e-7;        //relative error

        gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);
        //forward integration is usually easy:
        //— Function: int gsl_odeiv2_driver_apply (gsl_odeiv2_driver * d, double * t, const double t1, double y[])
        //  This function evolves the driver system d from t to t1. Initially vector y should contain the values of dependent variables at point t. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply is returned, and t and y contain the values from last successful step.
        
        //If maximum number of steps is reached, a value of GSL_EMAXITER is returned. If the step size drops below minimum value, the function returns with GSL_ENOPROG. If the user-supplied functions defined in the system sys returns GSL_EBADFUNC, the function returns immediately with the same return code. In this case the user must call gsl_odeiv2_driver_reset before calling this function again.
        
        int status = gsl_odeiv2_driver_apply (d, &t, tf, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
        }
        gsl_odeiv2_driver_free (d);
        
    }else{
        //backward integration.
        //sometimes, a nonstiff system might become stiff when integrating backward.
        //using a fixed time step.
        bool isNan=true;
        int order=0;
        double dist=1;
        while(isNan){// || (order<5 && dist>0.5)){
            for(int i=0;i<MAXEQ;i++) y[i]=state[i];
            double hstart=-1e-10;   //initial step size
            double epsabs=1e-7;     //absolute error
            double epsrel=1e-7;        //relative error
            gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);
            //— Function: int gsl_odeiv2_driver_apply_fixed_step (gsl_odeiv2_driver * d, double * t, const double h, const unsigned long int n, double y[])
            //This function evolves the driver system d from t with n steps of size h. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply_fixed_step is returned, and t and y contain the values from last successful step.
            //t===>0   -0 -1e-11
            cout << "t===>" << tf <<  "   " << -tf << " " << t0 << endl;
            int status =  gsl_odeiv2_driver_apply_fixed_step (d,&t, -tf, -1e-6, y);
            isNan=false;
            for(int i=0;i<MAXEQ;i++){
                if (isnan(y[i])) isNan=true ;
            }
            
            if(isNan) cout << "reducing integration time" << endl ;
            tf = tf/10.0;
            if (status != GSL_SUCCESS){
                printf ("error, return value=%d\n", status);
                isNan = true;
            }
            gsl_odeiv2_driver_free (d);
            order++;
            
            dist=0;
            for(int i=0;i<MAXEQ;i++){
                cout << y[i] << " " << state[i] << " " << y[i]-state[i] << endl ;
                dist += (y[i]-state[i])*(y[i]-state[i]);
            }
            cout << "dist="<< dist << endl ;
            dist=sqrt(dist);
            
        }
        
    }
    
    double* result = new double[MAXEQ];
    for(int i=0;i<MAXEQ;i++){
        result[i] = y[i];
    }
    return result ;
}

double* Circuit::integration_stiff(gsl_odeiv2_system sys, double* state, Configuration* config, bool direction, double t0){
    const gsl_odeiv2_step_type *steppers = gsl_odeiv2_step_bsimp;//gsl_odeiv2_step_bsimp; //gsl_odeiv2_step_bsimp;
    // initial values
    double y[MAXEQ];
    for(int i=0;i<MAXEQ;i++) y[i]=state[i];
    
    double dist=1; int order=0;
    double start = t0;
    double initstepsize = (direction)?1e-10:1e-10;
    double max_step = 1e-6;
    if(direction){
        config->getParameter("edu.uiuc.crhc.core.simulation.step", &max_step);
    }else{
        config->getParameter("edu.uiuc.crhc.core.simulation.backwardstep", &max_step);
    }
    max_step = 1e-6;
    
    start = t0;
    double end = start + ( (direction)?max_step:-max_step ) ;
    const double epsabs = 1e-15 ;
    const double epsrel = 1e-15 ;
    
    int s = 0;
    int steps = 0;
    double h = (direction)?initstepsize:-initstepsize;
    
    for(int i=0;i<MAXEQ;i++) y[i]=state[i];
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, h, epsabs, epsrel);
    int status = gsl_odeiv2_driver_apply (d, &start, end, y);
    if (status != GSL_SUCCESS){
        printf ("error, return value=%d\n", status);
    }

    
    /*
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_standard_new (&sys, steppers, h, epsabs, epsrel, 1.0, 0.0);
    gsl_odeiv2_driver_set_hmax(d, 1e-9);
    while ( direction?(start < end):(start>end)){
        s = gsl_odeiv2_evolve_apply (d->e, d->c, d->s, &sys, &start, end, &h, y);
        if (s != GSL_SUCCESS){
            gsl_test (s, "sys_driver: %s evolve_apply returned %d", gsl_odeiv2_step_name (d->s), s);
            break;
        }
        if (steps > 1e7){
            gsl_test (GSL_EMAXITER, "sys_driver: %s evolve_apply reached maxiter at t=%g", gsl_odeiv2_step_name (d->s), start);
            s = GSL_EMAXITER;
            break;
        }
        steps++;
    }*/
        
    gsl_odeiv2_step_reset(d->s);
    gsl_odeiv2_driver_free (d);
     
    
    double* result = new double[MAXEQ];
    for(int i=0;i<MAXEQ;i++){
        result[i] = y[i];
    }
    
    return result ;
}

//This function is the interface to the simulation, nothing special here.
// it just constructs the system based on the rhs-odes and select the solver based on stiff/nonstiff system
double* Circuit::integration(double* state, Configuration* config, bool direction, double t0){
    vector<double> p;
    double var_min = -0.1 ;
    double var_max = +0.1 ;
    for(int i=0;i<MAXEQ;i++){
        p.push_back( random(var_min, var_max) ) ;
    }
    
    gsl_odeiv2_system sys;
    switch (type) {
        case vanderpol:
            sys = {vanderpol_func, vanderpol_jac, 2, &p};
            return integration_nonstiff( sys, state, config, direction, t0);
        case tunneldiode:
            sys = {tunneldiode_func, tunneldiode_jac, 2, &p};
            return integration_nonstiff( sys, state, config, direction, t0);
        case ringmodulator:
            sys = {rhs_ringmod, jac_ringmod, 15, &p};
           //return integration_nonstiff( sys, state, config, direction, t0);

            return integration_stiff( sys, state, config, direction, t0);
            
        case dae_amplifier:
            return integration_dae(state, config, direction, t0);
        default:
            cout << "undefined circuit" << endl ;
            break;
    }
    
}

Trace* Circuit::simulate_stiff (bool direction){
    Trace* trace = new Trace();
    vector<double> p;
    for(int i=0;i<MAXEQ;i++){
        p.push_back( 0 ) ;
    }
    
    const gsl_odeiv2_step_type *steppers;
    
    // Required error tolerance for each stepper.
    double err_target;
    
    // initial values for each ode-solver
    double y[MAXEQ];
    
    
    // Problems, their names and number of equations
    const gsl_odeiv2_system sys = {rhs_ringmod, jac_ringmod, MAXEQ, &p};
    const char *probname = "ringmod" ;
    
    // Integration interval for problems
    const double start = 0;
    const double end = direction?1e-6:-1e-6 ;
    
    const double epsabs = 1e-7 ;
    const double epsrel = 1e-7 ;
    const double initstepsize = direction?1e-10:-1e-10 ;
    
    // Steppers
    //steppers = gsl_odeiv2_step_bsimp;
    err_target = 1e-3;
    steppers = gsl_odeiv2_step_msbdf;
    //err_target = 1e-3;
    
    //set the initial value
    for (int i = 0; i < MAXEQ; i++){
        y[i] = 0.0;
    }
    
    // This function evolves a system sys with stepper T from t0 to t1.
    // Step length is varied via error control with possibly different
    // absolute and relative error tolerances.
    int s = 0;
    int steps = 0;
        
    double t = start;
    double h = initstepsize;
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_standard_new (&sys, steppers, h, epsabs, epsrel, 1.0, 0.0);
        
    extern int nfe, nje;
    nfe = 0;
    nje = 0;
        
    cout <<"dir=" << direction << " " << t << " " << end << endl ;
    while (direction?(t < end):(t>end)){
        cout << steps << " " << t << endl ;
        cout << y[0] << " " << y[1] << endl ; 
        s = gsl_odeiv2_evolve_apply (d->e, d->c, d->s, &sys, &t, end, &h, y);
        if (s != GSL_SUCCESS){
            gsl_test (s, "sys_driver: %s evolve_apply returned %d", gsl_odeiv2_step_name (d->s), s);
            break;
        }
        if (steps > 1e7){
            gsl_test (GSL_EMAXITER, "sys_driver: %s evolve_apply reached maxiter at t=%g", gsl_odeiv2_step_name (d->s), t);
            s = GSL_EMAXITER;
            break;
        }
        steps++;
        vector<double> v;
        for(int i=0;i<MAXEQ;i++) v.push_back(y[i]) ;
        trace->AddSample(v);
    }
    gsl_test (s, "%s %s [%g,%g], %d steps (nfe %d, nje %d) completed", gsl_odeiv2_step_name (d->s), probname, start, end, steps, nfe, nje);
    gsl_odeiv2_driver_free (d);
        
    if (s != GSL_SUCCESS){
        printf ("start=%.5e, initstepsize=%.5e\n", start, initstepsize);
        gsl_test (s, "test_extreme_problems %s %s", steppers->name, probname);
    }
    
    return trace;
}

Trace* Circuit::simulate(Configuration* config, bool direction){
    Trace* trace = new Trace();

    int steps = 10000 ;
    vector<double> p;
    for(int i=0;i<MAXEQ;i++) p.push_back(0);
    
    gsl_odeiv2_system sys;
    switch (type) {
        case vanderpol:
            sys = {vanderpol_func, vanderpol_jac, 2, &p};
            break;
        case tunneldiode:
            sys = {tunneldiode_func, tunneldiode_jac, 2, &p};
            break;
        case ringmodulator:
            sys = {rhs_ringmod, jac_ringmod, 15, &p};
            break;
        default:
            cout << "undefined circuit" << endl ;
            break;
    }
    
    double hstart=direction?1e-10:-1e-10;
    double epsabs=0;
    double epsrel=1e-7;
    
    //gsl_odeiv2_step_rk2   Explicit embedded Runge-Kutta (2, 3) method.
    //gsl_odeiv2_step_rk4   Explicit 4th order (classical) Runge-Kutta. Error estimation is carried out by the step doubling method.
    //gsl_odeiv2_step_rkf45 Explicit embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a good general-purpose integrator.
    //gsl_odeiv2_step_rkck  Explicit embedded Runge-Kutta Cash-Karp (4, 5) method.
    //gsl_odeiv2_step_rk8pd Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.
    //gsl_odeiv2_step_rk1imp  Implicit Gaussian first order Runge-Kutta (implicit Euler or backward Euler method). Error estimation is carried out by the step doubling method. This algorithm requires the Jacobian and access to the driver object via gsl_odeiv2_step_set_driver.
    //gsl_odeiv2_step_rk2imp  Implicit Gaussian second order Runge-Kutta. Also known as implicit mid-point rule. Error estimation is carried out by the step doubling method. This stepper requires the Jacobian and access to the driver object via gsl_odeiv2_step_set_driver.
    //gsl_odeiv2_step_rk4imp    Implicit Gaussian 4th order Runge-Kutta. Error estimation is carried out by the step doubling method. This algorithm requires the Jacobian and access to the driver object via gsl_odeiv2_step_set_driver.
    //gsl_odeiv2_step_bsimp    Implicit Bulirsch-Stoer method of Bader and Deuflhard. The method is generally suitable for stiff problems. This stepper requires the Jacobian.
    //gsl_odeiv2_step_msadams        A variable-coefficient linear multistep Adams method in Nordsieck form. This stepper uses explicit Adams-Bashforth (predictor) and implicit Adams-Moulton (corrector) methods in P(EC)^m functional iteration mode. Method order varies dynamically between 1 and 12. This stepper requires the access to the driver object via gsl_odeiv2_step_set_driver.
    //gsl_odeiv2_step_msbdf A variable-coefficient linear multistep backward differentiation formula (BDF) method in Nordsieck form. This stepper uses the explicit BDF formula as predictor and implicit BDF formula as corrector. A modified Newton iteration method is used to solve the system of non-linear equations. Method order varies dynamically between 1 and 5. The method is generally suitable for stiff problems. This stepper requires the Jacobian and the access to the driver object via gsl_odeiv2_step_set_driver.

    const gsl_odeiv2_step_type *stepper = gsl_odeiv2_step_rk4imp;
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, stepper, hstart, epsabs, epsrel);

    double t = 0, t1 = direction?1e-3:-1e-3;
    double* state = getInitialState();
    double y[MAXEQ] ;
    for(int i=0;i<MAXEQ;i++) y[i]=state[i];

    y[0]=0.266439;
    y[1]=0.198983;
    y[2]=0.476002;
    y[3]=-0.182432;
    y[4]=-0.211285;
    y[5]=0.447149;
    y[6]=0.11102;
    y[7]=-2.86006e-06;
    y[8]=-1.30403e-06;
    y[9]=0.0011262;
    y[10]=0.000403679;
    y[11]=-0.00112524;
    y[12]=-0.000404643;
    y[13]=0.000386292;
    y[14]=-0.000308394;
    
    vector<double> v0 ;
    for(int i=0;i<dim;i++) v0.push_back(y[i]);
    trace->AddSample(v0);
    
    for (int i = 1; i <= steps; i++){
        cout << i << endl ;
        double ti = i * t1 / (double)steps;
        //int gsl_odeiv2_driver_apply_fixed_step (gsl_odeiv2_driver * d, double * t, const double h, const unsigned long int n, double y[])     This function evolves the driver system d from t with n steps of size h. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply_fixed_step is returned, and t and y contain the values from last successful step.
        int h= 1000;
        double step_size = -ti/h;
        //int status = gsl_odeiv2_driver_apply_fixed_step(d, &t, step_size, h, y);
        cout << t << " " << ti << endl;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        vector<double> v;
        for(int i=0;i<dim;i++) v.push_back(y[i]) ;
        trace->AddSample(v);
        cout << v[0] << " " << v[1] << " " ;
    }
    for(int i=0;i<dim;i++){ cout << y[i] << " " ; }
    cout << endl ;
    gsl_odeiv2_driver_free (d);
    return trace;
}

Trace* Circuit::backward_simulate(Configuration* config){
    Trace* trace = new Trace();
    int steps = 100000 ;
    vector<double> p;
    p.push_back(0);
    p.push_back(0);
    
    gsl_odeiv2_system sys;
    switch (type) {
        case vanderpol:
            sys = {vanderpol_func, vanderpol_jac, 2, &p};
            break;
        case tunneldiode:
            sys = {tunneldiode_func, tunneldiode_jac, 2, &p};
            break;
        case ringmodulator:
            sys = {rhs_ringmod, jac_ringmod, 15, &p};
            break;
        default:
            cout << "undefined circuit" << endl ;
            break;
    }
    
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double t = 0.0, t1 = 1e-7;
    cout << t1 << endl ;
    double* state = getInitialState();
    double y[15] = {state[0],
        state[1],
        state[2],
        state[3],
        state[4],
        state[5],
        state[6],
        state[7],
        state[8],
        state[9],
        state[10],
        state[11],
        state[12],
        state[13],
        state[14]};
    
    for(int i=0;i<dim; i++) y[i] = state[i];
    //    double y[2] = { state[0], state[1] };
    vector<double> v0 ;
    for(int i=0;i<dim;i++) v0.push_back(y[i]);
    trace->AddSample(v0);
    
    for (int i = 1; i <= steps; i++){
        double ti = i * t1 / (double)steps;
        /*
        — Function: int gsl_odeiv2_driver_apply (gsl_odeiv2_driver * d, double * t, const double t1, double y[])
        This function evolves the driver system d from t to t1. Initially vector y should contain the values of dependent variables at point t. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply is returned, and t and y contain the values from last successful step.
        
        If maximum number of steps is reached, a value of GSL_EMAXITER is returned. If the step size drops below minimum value, the function returns with GSL_ENOPROG. If the user-supplied functions defined in the system sys returns GSL_EBADFUNC, the function returns immediately with the same return code. In this case the user must call gsl_odeiv2_driver_reset before calling this function again.
        
        — Function: int gsl_odeiv2_driver_apply_fixed_step (gsl_odeiv2_driver * d, double * t, const double h, const unsigned long int n, double y[])
        This function evolves the driver system d from t with n steps of size h. If the function is unable to complete the calculation, an error code from gsl_odeiv2_evolve_apply_fixed_step is returned, and t and y contain the values from last successful step.
        */
        
        //int status = gsl_odeiv2_driver_apply (d, &t, -ti, y);
        int status =  gsl_odeiv2_driver_apply_fixed_step (d,&t, -ti, 100, y);
        cout << y[0] << " " << y[1] << endl ;
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        vector<double> v;
        for(int i=0;i<dim;i++){
            v.push_back(y[i]);
        }
        trace->AddSample(v);
    }
    gsl_odeiv2_driver_free (d);
    return trace;
}


Trace* Circuit::sim_test(){
    Trace* trace = new Trace();
    vector<double> p;
    for(int i=0;i<MAXEQ;i++){
        p.push_back( 0 ) ;
    }
    
    const gsl_odeiv2_step_type *stepper = gsl_odeiv2_step_msbdf;//gsl_odeiv2_step_bsimp ; //gsl_odeiv2_step_msbdf;
    
    //set the initial value
    double y[MAXEQ];
    for (int i = 0; i < MAXEQ; i++){
        y[i] = 0.0;
    }
    
    // Problems, their names and number of equations
    gsl_odeiv2_system sys = {rhs_ringmod, jac_ringmod, MAXEQ, &p};
 
    // Integration interval for problems
    double start = 0;
    double end = 1e-4 ;
    
    double epsabs = 1e-9 ;
    double epsrel = 1e-9 ;
    double initstepsize = 1e-10 ;
    double steps = 10000;
    
    cout << "-----------------Forward Sim--------------------------------" << endl ;
    
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkck , initstepsize, epsabs, epsrel);
    
    vector<double> v0 ;
    for(int i=0;i<dim;i++) v0.push_back(y[i]);
    trace->AddSample(v0);
    
    double t=0;
    gsl_odeiv2_driver_set_hmax(d, 1e-10);
    for (int i = 1; i <= steps; i++){
        double ti = i * end / (double)steps;
        double h_fixed_step_size=1e-11;
        if( h_fixed_step_size > end/steps ) cout << "step size mismatch" << endl ;
        double n_number_of_steps =(end / (double)steps) /  h_fixed_step_size   ;
        cout << "number of steps=" << n_number_of_steps << endl ;
        double n = 1 ;
        cout << i << " " << ti << endl ;
        //int status = gsl_odeiv2_driver_apply (d, &t, -ti, y);
        int status =  gsl_odeiv2_driver_apply_fixed_step (d, &t, h_fixed_step_size, n_number_of_steps, y);
        cout << t << "   ->" << y[0] << " " << y[1] << endl ;
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        vector<double> v;
        for(int i=0;i<dim;i++){
            v.push_back(y[i]);
        }
        trace->AddSample(v);
    }
    gsl_odeiv2_driver_free (d);
    
     cout << "-----------------Backward Sim--------------------------------" << endl ;
    // Integration interval for problems
    start = 0;
    end = 1e-4 ;
    
    epsabs = 1e-9 ;
    epsrel = 1e-9 ;
    initstepsize = -1e-10 ;
    steps = 10000;
    
    d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkck , initstepsize, epsabs, epsrel);
    
    t=1e-4;
    //gsl_odeiv2_driver_set_hmax(d, 1e-10);
    for (int i = 1; i <= steps; i++){
        double ti = i * end / (double)steps;
        double h_fixed_step_size=1e-11;
        if( h_fixed_step_size > end/steps ) cout << "step size mismatch" << endl ;
        double n_number_of_steps =( end / (double)steps) /  h_fixed_step_size   ;
        cout << "number of steps=" << n_number_of_steps << endl ;
        cout << i << " " << ti << endl ;
        //int status = gsl_odeiv2_driver_apply (d, &t, -ti, y);
        int status =  gsl_odeiv2_driver_apply_fixed_step (d, &t, -h_fixed_step_size, n_number_of_steps, y);
        cout << t << "   ->" << y[0] << " " << y[1] << endl ;
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        vector<double> v;
        for(int i=0;i<dim;i++){
            v.push_back(y[i]);
        }
        trace->AddSample(v);
    }
    gsl_odeiv2_driver_free (d);
    
    return trace;
}


Trace* Circuit::simulate_dae(){
    Trace* trace = new Trace();
    
    // dimension of problem
	const int n(8);
	// initial values for y
	double y[n] = {0.0, 6.0, 3.0, 3.0, 6.0, 3.0, 3.0, 0.0};
	//double y[n] = {0.24256355,   5.7655593 ,  2.3201511 ,  2.2570614,   3.8214509 ,  2.8499236 ,  3.0064879 , -0.0055545641};
	//double y[n] = {0.24841494,   5.7604951 ,  2.4532016 , -12057.931,  -12041.745 , -30.242248  ,-29.846618  ,-33.390254 };
	// initial value for x
	double xbeg(0);	//double xbeg(0.0);
	// final value for x
	const double xend(1); //const double xend(0.03);
	// interval of x for printing output
	double dx(0.00025);//double dx(0.0025);
	// rtoler and atoler are scalars
	int itoler(0);
	// relative tolerance
	//double *rtoler = new double(1.0e-5);
	double *rtoler = new double(1.0e-5);
	// absolute tolerance
	//double *atoler = new double(1.0e-11);
	double *atoler = new double(1.0e-11);
	// use SolutionOutput routine
	const int iout(1);
	// analytical Jacobian function provided
	const int ijac(1);
	// number of non-zero rows below main diagonal of Jacobian
	int mljac(1);
	// number of non-zero rows above main diagonal of Jacobian
	int mujac(2);
	// Mass matrix routine is provided
	const int imas(1);
	// number of non-zero rows below main diagonal of mass matrix
	int mlmas(1);
	// number of non-zero rows above main diagonal of mass matrix
	int mumas(1);
    
    // Use default values (see header files) for these parameters:
	double hinit(0.0);
	double hmax(0.0);
	int nmax(0);
	double uround(0.0), safe(0.0), facl(0.0), facr(0.0);
	int nit(0);
	bool startn(false);
	int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
	bool hess(false);
	double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.0);
    
    double steps=10000;
    for(int i=0;i<steps;i++){
        variation.clear();
        variation.push_back( random(-0.00001, 0.00001));
        double t0 = i*(xend-xbeg)/steps;
        double t1 = t0 + (xend-xbeg)/steps ;
        cout << "int from " << t0 << " to " << t1 << endl ;
        StiffIntegratorT stiffT(n, y, t0, t1, dx, itoler, rtoler, atoler,
                                iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
                                mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
                                m1, m2, hess, fnewt, quot1, quot2, thet);
        
        cout << "\n\n*******Problem integrated with RADAU5*******\n\n";
        
        stiffT.Integrate();
        
        // print statistics
        cout << "fcn = " << stiffT.NumFunction() <<
        " jac = " << stiffT.NumJacobian() <<
        " step = " << stiffT.NumStep() <<
        " accpt = " << stiffT.NumAccept() <<
        " rejct = " << stiffT.NumReject() <<
        " dec = " << stiffT.NumDecomp() <<
        " sol = " << stiffT.NumSol() << endl;

        
        vector<double> v;

        for(int j=0;j<8;j++){
            cout << y[j] << " . " ;
            v.push_back(y[j]);
        }
        trace->AddSample(v);

    }
    return trace;
}


double* Circuit::integration_dae(double* state, Configuration* config, bool direction, double t0){
    // initial values
    double y[8];
    for(int i=0;i<8;i++) y[i]=state[i];
    
    // dimension of problem
	const int n(8);
	// initial value for x
	double xbeg(t0);	//double xbeg(0.0);
    const double dt=1e-6;
	// final value for x
	const double xend(  xbeg + dt ); //const double xend(0.03);
	// interval of x for printing output
	double dx(0.00001);//double dx(0.0025);
	// rtoler and atoler are scalars
	int itoler(0);
	// relative tolerance
	//double *rtoler = new double(1.0e-5);
	double *rtoler = new double(1.0e-7);
	// absolute tolerance
	//double *atoler = new double(1.0e-11);
	double *atoler = new double(1.0e-14);
	// use SolutionOutput routine
	const int iout(1);
	// analytical Jacobian function provided
	const int ijac(1);
	// number of non-zero rows below main diagonal of Jacobian
	int mljac(1);
	// number of non-zero rows above main diagonal of Jacobian
	int mujac(2);
	// Mass matrix routine is provided
	const int imas(1);
	// number of non-zero rows below main diagonal of mass matrix
	int mlmas(1);
	// number of non-zero rows above main diagonal of mass matrix
	int mumas(1);
    
    // Use default values (see header files) for these parameters:
	double hinit(0.0);
	double hmax(0.0);
	int nmax(0);
	double uround(0.0), safe(0.0), facl(0.0), facr(0.0);
	int nit(0);
	bool startn(false);
	int nind1(0), nind2(0), nind3(0), npred(0), m1(0), m2(0);
	bool hess(false);
	double fnewt(0.0), quot1(0.0), quot2(0.0), thet(0.0);
    
        variation.clear();
        variation.push_back( random(-0.00001, 0.00001));
        StiffIntegratorT stiffT(n, y, xbeg, xend, dx, itoler, rtoler, atoler,
                                iout, hinit, hmax, nmax, uround, safe, facl, facr, ijac, mljac,
                                mujac, imas, mlmas, mumas, nit, startn, nind1, nind2, nind3, npred,
                                m1, m2, hess, fnewt, quot1, quot2, thet);
        
        cout << "\n\n*******Problem integrated with RADAU5*******\n\n";
        
        int status = stiffT.Integrate();
        
        // print statistics
        cout << "fcn = " << stiffT.NumFunction() <<
        " jac = " << stiffT.NumJacobian() <<
        " step = " << stiffT.NumStep() <<
        " accpt = " << stiffT.NumAccept() <<
        " rejct = " << stiffT.NumReject() <<
        " dec = " << stiffT.NumDecomp() <<
        " sol = " << stiffT.NumSol() << endl;
        
        
    
    double* result = new double[8];
    for(int i=0;i<8;i++){
        result[i] = y[i];
    }
    
    return result ;
}





