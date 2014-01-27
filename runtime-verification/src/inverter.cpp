#include "inverter.h"
#include <fstream>
#include <iostream>
#include <sstream>

vector<double>  Inverter::simulate(double* ic, vector<double> param, vector<string> setting, double dt){
	vector<double> result;
	
	stringstream sed; 
	sed << "cat inverter_template.sp | sed ";
	for (int i = 0; i < param.size(); i++){
		sed << "-e s/$PARAM_" << i << "/" << param[i] << "/ ";
	}
	sed << "-e s/$DT/" << dt << "/ ";
	sed << "-e s/$SAVE_FILE/" << setting[0] << "/ ";
	sed << "-e s/$LOAD_FILE/" << setting[1] << "/ ";
	sed << " > inverter_netlist.sp" << endl;

	cout << sed.str() << endl;
	system(sed.str().c_str());
	system("hspice inverter_netlist.sp > Sim.txt");
	system("cat Sim.txt | grep 5.0000000000n > grep.txt");		// grep the line containing my results at sim time 10ns

	string line;
	ifstream simResult("grep.txt");
	while (simResult.good()){
		getline(simResult, line);
		if (line.size() > 0){
			cout << "Sim: " << line << endl;
			result = System::parse(line);
			break;
		}
	}
	simResult.close();
	cout << "Here! " << endl;
	cout << result.size() << endl;
	return result;
}