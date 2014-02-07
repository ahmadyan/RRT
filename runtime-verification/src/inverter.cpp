#include "inverter.h"
#include <fstream>
#include <iostream>
#include <sstream>

Inverter::Inverter(Configuration* config){
	cout << "Selected system is Inverter" << endl;
	config->getParameter("edu.uiuc.csl.system.dimension", &d);
	System::simulator = SPICE;
}

vector<double>  Inverter::simulate(double* ic, vector<double> param, vector<string> setting, double t0, double dt){
	vector<double> result;
	
	bool dcSimulation = false;
	if (setting[2].compare("dc") == 0){
		dcSimulation = true;
	}

	stringstream sed; 
	sed << "cat inverter_template.sp | sed ";
	for (int i = 0; i < param.size(); i++){
		sed << "-e s/$PARAM_" << i << "/" << param[i] << "/ ";
	}
	if (dcSimulation){		
		sed << "-e s/$tran/" << "\".tran 1fs " << dt << "\"/ ";
		sed << "-e s/$save/" << "\".save type=ic file=" << setting[0] << " level=all time=" << dt << "\"/ ";
	}else{
		sed << "-e s/$tran/" << "\".tran 1fs " << dt << " uic" << "\"/ ";
		sed << "-e s/$save/" << "\".save type=ic file=" << setting[0] << " level=all time=" << dt << "\"/ ";
		sed << "-e s/$load/" << "\".load file=" << setting[1] << "\"/";
	}
	sed << " > inverter_netlist.sp" << endl;

	cout << sed.str() << endl;
	system(sed.str().c_str());
	system("hspice inverter_netlist.sp > Sim.txt");

	//system("cat Sim.txt | grep 5.0000000000n > grep.txt");		// grep the line containing my results at sim time 10ns
	//string line;
	//ifstream simResult("grep.txt");
	//while (simResult.good()){
	//	getline(simResult, line);
	//	if (line.size() > 0){
	//		cout << "Sim: " << line << endl;
	//		result = System::parse(line);
	//		break;
	//	}
	//}
	//simResult.close();
	result.push_back(dt);
	double* tmp = parseICFile(setting[0] + "0");
	for (int i = 0; i < getDimension(); i++){
		result.push_back(tmp[i]);
	}
	return result;
}