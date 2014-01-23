#include "System.h"
#define USE_HSPICE	true
System::System(CircuitType t){
    type=t;
}

System::System(){}
System::~System(){}

void System::setSystem(int _d,int (*f)(double t, const double y[], double f[], void *params), 
                    int (*j) (double t, const double y[], double *dfdy, double dfdt[], void *params) ){
    d=_d;
    function = f ;
    jacobian = j ;
}

void System::dumpSystemInfo(){
}

//this is a simple integrator that would do a step integration
double System::simulate(double* initialState, double* param, double dt){
	if(USE_HSPICE)
		return runSPICE(initialState, param, dt);
	else
		return integareODE(initialState, param, dt);
}

double System::integareODE(double* initialState, double* param, double dt){
/*     double mu = param;
    gsl_odeiv2_system dynamic = {function, jacobian, d , &mu};  // {vector-element, vector-of-derivative, system-dimension, arbitrary-parameters}
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&dynamic, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double t=0;
    int status = gsl_odeiv2_driver_apply (d, &t, dt, initialState);
    if (status != GSL_SUCCESS){
        cout << "error, return value=" <<  status << endl ;
        return 0;
    }
    cout <<  t  << " " << initialState[0] << " " <<  initialState[1] << endl ;
    gsl_odeiv2_driver_free (d);
    return t;*/
	return 0 ;
}

double System::runSPICE(double* initialState, double* param, double dt){
	spice hspice;
	vector<double> result = hspice.simulate(type, initialState, param, dt);
	cout << "result.size()=" << result.size() << endl;
	for (int i = 0; i < result.size(); i++){
		cout << i << " " << result[i] << endl;
	}
	for(int i=0;i<result.size()-1;i++){
		initialState[i]=result[i+1];
		cout << i << " " << result[i+1];
	}
	return result[0];
}


//this is a compelete integrator, does not have any particular use, except for debugging the previous function
double System::simulate(){
/*    double mu = 10;
    gsl_odeiv2_system dynamic = {function, jacobian, 2, &mu};
    gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&dynamic, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);
    double t = 0.0, t1 = 100.0;
    double y[2] = { 1.0, 0.0 };
    for (int i = 1; i <= 100; i++){
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);
        if (status != GSL_SUCCESS){
            printf ("error, return value=%d\n", status);
            break;
        }
        printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }
    gsl_odeiv2_driver_free (d);
    return t1 ;*/
	return 0 ;
}

void System::setDimension(int _d){
    d=_d;
}

int System::getDimension(){
    return d;
}
