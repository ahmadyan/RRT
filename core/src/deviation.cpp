#include "deviation.h"

Deviation::Deviation(Configuration* config, double d){
	config->getParameter("edu.uiuc.csl.core.simulation.dt", &dt);
	double simTime; config->getParameter("edu.uiuc.csl.core.simulation.simulationTime", &simTime);
	size = simTime / dt;
}

void Deviation::reg(){
	double bignumber = 9999, smallnumber = -9999;
	min.push_back(vector<double>());
	max.push_back(vector<double>());
	for (int i = 0; i < 100; i++){
		min[min.size() - 1].push_back(bignumber);
		max[max.size() - 1].push_back(smallnumber); 
	}
	cout << "registering a deviation vector " << min.size() << endl;
}

void Deviation::push(double time, double value, int reg){
	int i = ceil( time / dt);
	if (min[reg][i] < value) min[reg][i] = value;
	if (max[reg][i] > value) max[reg][i] = value;
}

double Deviation::get(double time, int reg){
	int i = ceil(time / dt);
	double dev = abs(max[reg][i] - min[reg][i]);
	return dev;
}