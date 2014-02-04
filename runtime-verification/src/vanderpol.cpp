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

//Func will define the dynamics of the system. This is Van der Pol Oscillator
int vanderpol_func(double t, const double y[], double f[], void *params){
	vector<double> p = *(vector<double>*)params;
	double p0 = p.at(0);    //variation in diode
	double p1 = p.at(1);    //variation in voltage source
	double mu = 2;
	f[0] = y[1] + p0;
	f[1] = -(y[0] + p1) - mu*y[1] * (y[0] * y[0] - 1);
	return GSL_SUCCESS;
}

int vanderpol_jac(double t, const double y[], double *dfdy, double dfdt[], void *params){
	vector<double> p = *(vector<double>*)params;
	double p0 = p.at(0);    //variation in diode
	double p1 = p.at(1);    //variation in voltage source
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
vector<double>  Vanderpol::simulate(double* ic, vector<double> param, vector<string> setting, double dt){
	double v0 = ic[0];
	double i0 = ic[1]; 
	double dvin = param[0];
	double dId = param[1];
	vector<double> result;
	FILE * netList;
	string aString;
	double rise, fall, delay;
	char nextString[80];
	netList = fopen("tdo.sp", "w");

	fprintf(netList, "TDO - TUNNEL DIODE OSCILLATOR	\n");
	fprintf(netList, "VIN	3	0	%f\n", 0.3 + dvin);
	//fprintf(netList, "R1	2	3	0.2\n");	//For oscillation result
	fprintf (netList, "R1	2	3	0.5\n");	//For No oscillation result
	fprintf(netList, "LS  2 	1 	1UH	 IC=%f \n", i0);
	fprintf(netList, "CS  1 	0 	1000PF\n");
	fprintf(netList, "G1 1 0 	POLY(1)	1 0 %f 0.6 -1.5 1 	\n", dId);
	fprintf(netList, ".IC V(1)=%f\n", v0);
	//fprintf (netList, ".TRAN 10NS 20NS 0 5NS UIC\n");
	fprintf(netList, ".TRAN 1NS 5NS UIC\n");
	//fprintf (netList, ".PLOT TRAN V(1) \n");
	fprintf(netList, ".PRINT V(1) I(LS) I(Vin)\n");
	fprintf(netList, ".OPT BRIEF numdgt=10\n");
	fprintf(netList, ".END\n");

	fclose(netList);
	system("hspice tdo.sp > Sim.txt");
	system("cat Sim.txt | grep 5.0000000000n > grep.txt");		// grep the line containing my results at sim time 10ns

	string line;
	ifstream simResult("grep.txt");
	while (simResult.good()){
		getline(simResult, line);
		if (line.size() > 0){
			cout << "Sim: " << line << endl;
			result = System::parse(line);
			break;
		}
	}
	simResult.close();
	cout << result.size() << endl;
	return result;

	
}

