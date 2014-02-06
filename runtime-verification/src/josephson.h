
#include <vector>
#include <string>
#include <gsl/gsl_odeiv2.h>
#include "system.h"
using namespace std;

class Josephson: public System {
	double mu;
public:
	Josephson();
	vector<double> simulate(double* ic, vector<double> parameters, vector<string> settings, double t0, double dt);
};