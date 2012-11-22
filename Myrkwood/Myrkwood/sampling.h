#pragma once

#include "config.h"
#include "point.h"

class sampling {
	//enum samplingType {uniform, gaussian, goalOriented};
	//int d;
	//double* min;
	//double* max;
	//samplingType type;
public:
	sampling(Configuration* config);
	Point* sample();
	virtual ~sampling();
};
