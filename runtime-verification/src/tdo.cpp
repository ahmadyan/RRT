#include "tdo.h"
#include <fstream>
#include <iostream>


/* Sample Netlist:
TDO - TUNNEL DIODE OSCILLATOR
VIN	3	0	300mv
R1	2	3	0.2
LS  2 	1 	1UH	 IC=71.0419253060m
CS  1 	0 	1000PF
G1 1 0 	POLY(1)	1 0 0 0.6 -1.5 1
.IC V(1)=222.8756004648m *I(LS)=48.3587298573m
.TRAN 5NS 50US 0 5NS UIC
*.SAVE TYPE=IC level=all time=10ns
.PLOT TRAN V(1)
.PRINT V(1) I(LS) I(Vin)
.OPT LIST NODE OPTS  numdgt=10 post
.END
*/
vector<double>  TDO::simulate(double* ic, vector<double> param, vector<string> setting, double dt){
	double v0 = ic[0];
	double i0 = ic[1]; 
	double dvin = param[0];
	double dId = param[1];
	vector<double> result;
	FILE * netList;
	string aString;
	double rise, fall, delay;
	char nextString[80];
	netList = fopen("tdo.sp", "w");

	fprintf(netList, "TDO - TUNNEL DIODE OSCILLATOR	\n");
	fprintf(netList, "VIN	3	0	%f\n", 0.3 + dvin);
	//fprintf(netList, "R1	2	3	0.2\n");	//For oscillation result
	fprintf (netList, "R1	2	3	0.5\n");	//For No oscillation result
	fprintf(netList, "LS  2 	1 	1UH	 IC=%f \n", i0);
	fprintf(netList, "CS  1 	0 	1000PF\n");
	fprintf(netList, "G1 1 0 	POLY(1)	1 0 %f 0.6 -1.5 1 	\n", dId);
	fprintf(netList, ".IC V(1)=%f\n", v0);
	//fprintf (netList, ".TRAN 10NS 20NS 0 5NS UIC\n");
	fprintf(netList, ".TRAN 1NS 5NS UIC\n");
	//fprintf (netList, ".PLOT TRAN V(1) \n");
	fprintf(netList, ".PRINT V(1) I(LS) I(Vin)\n");
	fprintf(netList, ".OPT BRIEF numdgt=10\n");
	fprintf(netList, ".END\n");

	fclose(netList);
	system("hspice tdo.sp > Sim.txt");
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
	cout << result.size() << endl;
	return result;

	/*
	ifstream simResult ("Sim.txt");

	while (simResult.good()){
	getline (simResult,line);
	cout << line << endl ;
	if (std::string::npos != line.find("transient analysis")){
	std::cout << "found! " << line << std::endl;
	//Now gettings transient results ...
	for(int i=0;i<9;i++){
	getline (simResult,line);
	cout << i << line << endl ;
	result = parse(line);
	if(i==5) return result ;	// this is the 5ns sim. the 5th line after transient analysis.
	for(int j=0;j<result.size();j++){
	cout.precision(10);
	cout << j << " " << result[j] << endl ;
	}

	}
	}
	}
	simResult.close();

	return result;
	*/
}