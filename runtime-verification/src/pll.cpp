#include <stdio.h>
#include <fstream>
#include "pll.h"


PLL::PLL(Configuration* config){
	config->getParameter("edu.uiuc.csl.system.dimension", &d);
	System::simulator = SPICE;
}

vector<double>  PLL::simulate(double* nodeset, vector<double> variation, vector<string> settings, double t0, double dt){
	vector<double> result;
	result.push_back(dt);
	stringstream sed;
	sed << "cat pll_template.sp | sed ";
	for (int i = 0; i < variation.size(); i++){
		sed << "-e s/$PARAM_" << i << "/" << variation[i] << "/ ";
	}
	sed << "-e s/$DT/" << dt << "/ ";
	sed << "-e s/$SAVE_FILE/" << settings[0] << "/ ";
	sed << "-e s/$LOAD_FILE/" << settings[1] << "/ ";
	sed << " > pll_netlist.sp" << endl;
	cout << sed.str() << endl;
	system(sed.str().c_str());
	generateICFile(nodeset, settings[1]);
	system("hspice pll_netlist.sp > Sim.txt");
	double* tmp = parseICFile(settings[0] + "0");
	for (int i = 0; i < getDimension(); i++){
		result.push_back(tmp[i]);
	}
	delete tmp;
	return result;
}

void  PLL::setInitialPLLState(double* state){
	/*
	* .nodeset
	+ e =  -3.5000
	+ eb =  -3.3000
	+ in =   0.
	+ inb =   0.
	+ mout =  -3.5000
	+ moutb =  -3.5000
	+ osc =-999.7550m
	+ oscb =  -1.0002
	+ out =  -3.5000
	+ outb =  -3.3000
	+ xvco.c =   1.0000
	+ xvco.s =  61.2612u
	+ xvco.s_clip =  81.6816u
	+ xpd.clip1 =   0.
	+ xpd.clip2 = 653.4528u
	+ xpd.n1 =   0.
	*/

	state[pll_e] = -3.5000;
	state[pll_eb] = -3.3000;
	state[pll_in] = 0;
	state[pll_inb] = 0;
	state[pll_mout] = -3.5000;
	state[pll_moutb] = -3.5000;
	state[pll_osc] = -999.7550e-3;
	state[pll_oscb] = -1.0002;
	state[pll_out] = -3.5000;
	state[pll_outb] = -3.3000;
	state[pll_xvco_c] = 1.0000;
	state[pll_xvco_s] = 61.2612e-6;
	state[pll_xvco_s_clip] = 81.6816e-6;
	state[pll_xpd_clip1] = 0;
	state[pll_xpd_clip2] = 653.4528e-6;
	state[pll_xpd_n1] = 0;
	state[pll_time] = 0;
}


void PLL::generateICFile(double* nodeset, string fileName){
	FILE *ic = fopen(fileName.c_str(), "w");

	fprintf(ic, ".option \n");
	fprintf(ic, "	+ gmindc=   1.0000p       \n");
	fprintf(ic, "	.nodeset	\n");
	fprintf(ic, "	+ e =  %f   \n", nodeset[pll_e]);
	fprintf(ic, "	+ eb =  %f       \n", nodeset[pll_eb]);
	fprintf(ic, "	+ in =   %f            \n", nodeset[pll_in]);
	fprintf(ic, "	+ inb =   %f            \n", nodeset[pll_inb]);
	fprintf(ic, "	+ mout =  %f        \n", nodeset[pll_mout]);
	fprintf(ic, "	+ moutb =  %f        \n", nodeset[pll_moutb]);
	fprintf(ic, "	+ osc = %f       \n", nodeset[pll_osc]);
	fprintf(ic, "	+ oscb =  %f        \n", nodeset[pll_oscb]);
	fprintf(ic, "	+ out =  %f        \n", nodeset[pll_out]);
	fprintf(ic, "	+ outb =  %f        \n", nodeset[pll_outb]);
	fprintf(ic, "	+ xvco.c =   %f        \n", nodeset[pll_xvco_c]);
	fprintf(ic, "	+ xvco.s =  %f       \n", nodeset[pll_xvco_s]);
	fprintf(ic, "	+ xvco.s_clip =  %f       \n", nodeset[pll_xvco_s_clip]);
	fprintf(ic, "	+ xpd.clip1 =   %f            \n", nodeset[pll_xpd_clip1]);
	fprintf(ic, "	+ xpd.clip2 = %f       \n", nodeset[pll_xpd_clip2]);
	fprintf(ic, "	+ xpd.n1 =   %f            \n", nodeset[pll_xpd_n1]);

	fclose(ic);
}