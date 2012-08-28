/*#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>

using namespace std;

double getDelay(double vTh1, double vTh2, double vTh3, double vTh4, double vTh5, double vTh6){

	FILE * netList;
	FILE * simResult;
	string aString;
	double rise, fall, delay;
	char nextString[80];
	//Line below can be commented out when you are passing these values into the function.
	//float vTh1=-0.0, vTh2=0.0, vTh3=-0.0, vTh4=0.0, vTh5=0.0, vTh6=0.0;
			
	netList = fopen("6TCell.cir", "w");

	fputs ("*6T SRAM cell\n.LIB '45nm_lib.pm' TT\n", netList);
	fputs ("MPL OUT IN VDD VDD pmos W=0.135000U L=0.050000U AS=0.009112P AD=0.009112P PS=0.270000U PD=0.270000U DELVTO=", netList);
	fprintf (netList, "%f\n", vTh1);
		
	fputs ("MNL OUT IN 0 0 nmos W=0.090000U L=0.050000U AS=0.004050P AD=0.004050P PS=0.180000U PD=0.180000U DELVTO=", netList);
	fprintf (netList, "%f\n", vTh2);
		
	fputs ("MPR IN OUT VDD VDD pmos W=0.135000U L=0.050000U AS=0.009112P AD=0.009112P PS=0.270000U PD=0.270000U DELVTO=", netList);
	fprintf (netList, "%f\n", vTh3);
	
	fputs ("MNR IN OUT 0 0 nmos W=0.090000U L=0.050000U AS=0.004050P AD=0.004050P PS=0.180000U PD=0.180000U DELVTO=", netList);
	fprintf (netList, "%f\n", vTh4);
		
	fputs ("MNL2 OUT WL BLN 0 nmos W=0.090000U L=0.050000U AS=0.004050P AD=0.004050P PS=0.180000U PD=0.180000U DELVTO=", netList);
	fprintf (netList, "%f\n", vTh5);
		
	fputs ("MNR2 IN WL BL 0 nmos W=0.090000U L=0.050000U AS=0.004050P AD=0.004050P PS=0.180000U PD=0.180000U DELVTO=", netList);
	fprintf (netList, "%f\n", vTh6);
		
	fputs ("VVDD VDD 0 1.1\n", netList);
	fputs ("VWL WL 0 PWL(0 1.1 7.5n 1.1 7.6n 0 15n 0 15.1n 1.1 19n 1.1 19.1n 0 24n 0 24.1n 1.1)\n", netList);
	fputs ("VBL BL 0 PULSE 0 1.1 0n 0 0 10n 20n\n", netList);
	fputs ("VBLN BLN 0 PULSE 0 1.1 10n 0 0 10n 20n\n", netList);
	fputs (".TRAN 1ps 35ns\n", netList);

	fputs (".MEASURE TRAN RISE_DELAY TRIG V(WL) VAL=0.55 TD=12n RISE=1 TARG V(OUT) VAL=0.55 TD=12n RISE=1\n", netList);
	fputs (".MEASURE TRAN FALL_DELAY TRIG V(WL) VAL=0.55 TD=23n RISE=1 TARG V(OUT) VAL=0.55 TD=23n FALL=1\n", netList);
	fputs (".MEASURE DELAY param='0.5*RISE_DELAY+0.5*FALL_DELAY'\n", netList);
	fputs (".END", netList);
	
        fclose(netList);
	//system("hspice 6t.cir");
        system ("hspice 6TCell.cir > Sim.txt");
	//System call to run script/simulation and store it in Sim.txt
        //system ("ps"); 	
	simResult = fopen("Sim.txt", "r");
	fseek (simResult, -1500, SEEK_END);
	
	do{
		fscanf (simResult, "%s", nextString);
		aString = "";
		for(int j = 0; nextString[j] != 0; j++)
			aString += nextString[j];
	} while(aString != "rise_delay=" && !feof(simResult));
	fscanf (simResult, "%lf", &rise);
	
	do{
		fscanf (simResult, "%s", nextString);
		aString = "";
		for(int j = 0; nextString[j] != 0; j++)
			aString += nextString[j];
	} while(aString != "fall_delay=" && !feof(simResult));
	fscanf (simResult, "%lf", &fall);

	do{
		fscanf (simResult, "%s", nextString);
		aString = "";
		for(int j = 0; nextString[j] != 0; j++)
			aString += nextString[j];
	} while(aString != "delay=" && !feof(simResult));
	fscanf (simResult, "%lf", &delay);

	//cout << "Value of delays: " << rise << " " << fall << " " << delay << endl;
	fclose(simResult);
        remove("Sim.txt");
  	//fclose(netList);
	//fclose(simResult);
    return delay;
}
*/