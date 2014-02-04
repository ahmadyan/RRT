
#include <vector>
#include <string>
#include "system.h"
using namespace std;

class Template: public System {
public:
	vector<double> simulate(double* ic, vector<double> parameters, vector<string> settings, double dt);
};