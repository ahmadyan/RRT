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

class Simulation{
public:
	virtual vector<double> simulate(vector<double> ic, vector<double> parameters, vector<string> settings, double dt);
	vector<double> parse(string str);
	bool is_only_ascii_whitespace(const std::string& str);
	double unit(char u);
};
