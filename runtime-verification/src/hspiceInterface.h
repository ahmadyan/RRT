#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <strstream>
#include <string>
#include <cctype>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
#include <vector>
using namespace std;

enum CircuitType {  TDO, //Tunnel Diode Oscillation
					PLL, //Phase-Locked loop circuit
					UNKNOWN };

class spice{
public:
	// I need --> Initial Condition ...
	// Simulation time
	// variation/input parameter
	// I return --> new state variable
	vector<double> simulate(CircuitType type, double* ic, double* param);
	vector<double> simulateTDO(double v0, double i0, double vin, double dId);
	vector<double> simulatePLL(double* ic, double param);
	vector<double> parse(string str);
};
