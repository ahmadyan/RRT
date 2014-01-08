#include "hspiceInterface.h"
#include "pll.h"
#define	unluckyThirteen	13

vector<double> spice::simulate(CircuitType type, double* ic, double* param){
	if(type==TDO){
		double v0 = ic[0];
		double i0 = ic[1];
		double dvin = param[0] ;
		double dId  = param[0]/10.0 ;
		return simulateTDO(v0,i0, dvin, dId);
	}else{
		return simulatePLL(ic, param[0]);
	}
	
	
}

bool is_only_ascii_whitespace( const std::string& str ){
	//bool isOnlyWhiteSpace = false;

	//return str.find_first_not_of (' ') == str.npos)
	std::string::const_iterator  it = str.begin();
    do {
        if (it == str.end()) return true;
    } while (*it >= 0 && *it <= 0x7f && std::isspace(*(it++)));
             // one of these conditions will be optimized away by the compiler,
             // which one depends on whether char is signed or not
    return false;
}

double unit(char u){
	switch(u){
	case 'f':
	case 'F':
		return 1e-15 ;
		break;
	case 'p':
	case 'P':
		return 1e-12 ;
		break;
	case 'n':
	case 'N':
		return 1e-9 ;
		break;
	case 'u':
	case 'U':
		return 1e-6 ;
		break;
	case 'm':
	case 'M':
		return 1e-3 ;
		break;
	case 'k':
	case 'K':
		return 1e3 ;
		break;

	default:
		return 1 ;
	}
}

//The input is something like this:
//5.0000000000n   226.6820055575m   71.3468319456m  -71.3468319456m
// The first number is transient simulation time, the rest are variables that need to be converted to double
vector<double> spice::parse(string s){
	vector<double> result ;
	string word ;
	string str = s ;
	if(s.c_str()[s.length()-1] == unluckyThirteen){	//guess what this does!
		str = s.substr(0, s.length()-1);
	}
	stringstream stream(str);
	while( getline(stream, word, ' ') ){
		if( !is_only_ascii_whitespace(word) ){
			double d ;
			stringstream ss(word);
			ss >> d ;
			char c = word.c_str()[word.length()-1] ;
			d = d*unit(c);
			result.push_back(d);
		}
	}
	return result ;
}


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
vector<double>  spice::simulateTDO(double v0, double i0, double dvin, double dId){
	vector<double> result;
	FILE * netList;
	string aString;
	double rise, fall, delay;
	char nextString[80];
	netList = fopen("tdo.sp", "w");

	fprintf (netList, "TDO - TUNNEL DIODE OSCILLATOR	\n");
	fprintf (netList, "VIN	3	0	%f\n", 0.3+dvin);
	fprintf (netList, "R1	2	3	0.2\n");	//For oscillation result
	//fprintf (netList, "R1	2	3	0.5\n");	//For No oscillation result
	fprintf (netList, "LS  2 	1 	1UH	 IC=%f \n", i0);
	fprintf (netList, "CS  1 	0 	1000PF\n");
	fprintf (netList, "G1 1 0 	POLY(1)	1 0 %f 0.6 -1.5 1 	\n", dId);
	fprintf (netList, ".IC V(1)=%f\n", v0);
	//fprintf (netList, ".TRAN 10NS 20NS 0 5NS UIC\n");
	fprintf (netList, ".TRAN 1NS 5NS UIC\n");
	//fprintf (netList, ".PLOT TRAN V(1) \n");
	fprintf (netList, ".PRINT V(1) I(LS) I(Vin)\n");
	fprintf (netList, ".OPT BRIEF numdgt=10\n");
	fprintf (netList, ".END\n");
	
	fclose(netList);
	system ("hspice tdo.sp > Sim.txt");
	system("cat Sim.txt | grep 5.0000000000n > grep.txt");		// grep the line containing my results at sim time 10ns

	string line;
	ifstream simResult("grep.txt");
	while (simResult.good()){
		getline(simResult, line);
		if (line.size() > 0){
			cout << "Sim: " << line << endl;
			result = parse(line);
			break;
		}
	}
	simResult.close();
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

vector<double>  spice::simulatePLL(double* nodeset, double variation){
	vector<double> result ;
	double dt = 5e-6;
	result.push_back(dt);
	FILE * ic;
	string aString;
	double rise, fall, delay;
	char nextString[80];
	ic = fopen("pll.ic0", "w");

	fprintf (ic, ".option \n");
	fprintf (ic, "	+ gmindc=   1.0000p       \n");
	fprintf (ic, "	.nodeset	\n");
	fprintf (ic, "	+ e =  %f   \n", nodeset[pll_e]);
	fprintf (ic, "	+ eb =  %f       \n", nodeset[pll_eb]);
	fprintf (ic, "	+ in =   %f            \n" , nodeset[pll_in]);
	fprintf (ic, "	+ inb =   %f            \n", nodeset[pll_inb]);
	fprintf (ic, "	+ mout =  %f        \n", nodeset[pll_mout]);
	fprintf (ic, "	+ moutb =  %f        \n", nodeset[pll_moutb]);
	fprintf (ic, "	+ osc = %f       \n", nodeset[pll_osc]);
	fprintf (ic, "	+ oscb =  %f        \n", nodeset[pll_oscb]);
	fprintf (ic, "	+ out =  %f        \n", nodeset[pll_out]);
	fprintf (ic, "	+ outb =  %f        \n", nodeset[pll_outb]);
	fprintf (ic, "	+ xvco.c =   %f        \n", nodeset[pll_xvco_c]);
	fprintf (ic, "	+ xvco.s =  %f       \n", nodeset[pll_xvco_s]);
	fprintf (ic, "	+ xvco.s_clip =  %f       \n", nodeset[pll_xvco_s_clip]);
	fprintf (ic, "	+ xpd.clip1 =   %f            \n", nodeset[pll_xpd_clip1]);
	fprintf (ic, "	+ xpd.clip2 = %f       \n", nodeset[pll_xpd_clip2]);
	fprintf (ic, "	+ xpd.n1 =   %f            \n", nodeset[pll_xpd_n1]);

	fclose(ic);
	ic = fopen("param.sp", "w");
	fprintf (ic, ".param var=%f\n", variation);
	fclose(ic);

	system ("hspice pll.sp > Sim.txt");

	string line;
	ifstream simResult ("pll.ic0");

	if (simResult.good()){
		for(int i=0;i<12;i++)
			getline (simResult,line);

		for(int i=0;i<16;i++){
			getline (simResult,line);
			line = line.substr( line.find_first_of("=")+1, line.length());
			double d ;
			stringstream ss(line);
			ss >> d ;
			line.erase(line.find_last_not_of(" \n\r\t")+1);
			char c = line.c_str()[line.length()-1] ;
			d = d*unit(c);
			result.push_back(d);
		}
    }
    simResult.close();

	return result;
}
