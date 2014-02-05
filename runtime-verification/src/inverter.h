
#include <vector>
#include <string>
#include "system.h"
using namespace std;

class Inverter: public System {
public:
	Inverter();
	vector<double> simulate(double* ic, vector<double> parameters, vector<string> settings, double t0, double dt);
};