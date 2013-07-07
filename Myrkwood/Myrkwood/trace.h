#pragma once
#include <vector>
using namespace std;
class Trace{
    vector< vector<double> > data;
public:
	Trace();
	~Trace();
    void AddSample(vector<double> sample);
    vector<double> getSample(int i);
    int getSampleSize();
    vector< vector<double> > getData();
};