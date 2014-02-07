
#include <vector>
#include <string>
#include "system.h"
using namespace std;

class TDO: public System {
public:
	TDO(Configuration* config);
	vector<double> simulate(double* ic, vector<double> parameters, vector<string> settings, double t0 , double dt);
};