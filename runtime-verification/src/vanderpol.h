
#include <vector>
#include <string>
#include <gsl/gsl_odeiv2.h>
#include "system.h"
using namespace std;

class Vanderpol: public System {
	double mu;
public:
	Vanderpol();
	vector<double> simulate(double* ic, vector<double> parameters, vector<string> settings, double dt);
};