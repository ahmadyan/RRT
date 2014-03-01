#pragma once
#include "System.h"
#include <iostream>
#include <string>
#include <vector>
using namespace std;

class OdeSolver : public System {
public:
    OdeSolver();
    ~OdeSolver();
	vector<double>  simulate(double* ic, vector<double> param, vector<string> setting, double dt);
};
