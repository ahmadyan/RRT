
#include <vector>
#include <string>
#include "system.h"
using namespace std;

class SPICE : public System {
public:
	SPICE(Configuration* config);
	vector<double> simulate(double* ic, vector<double> parameters, vector<string> settings, double t0, double dt);
};