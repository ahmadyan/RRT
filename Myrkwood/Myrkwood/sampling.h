#pragma once

#include "circuit.h"
#include "point.h"

//Random Point Generator
class RPG {
	enum samplingType {uniform, gaussian, goalOriented};
    Circuit* circuit;
public:
	RPG(Circuit* config);
    geometry::Point* sample();
	~RPG();
    double random(double a, double b);
};
