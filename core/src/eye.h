#pragma once

#include <iostream>
#include <vector>
#include "node.h" 

using namespace std;

class EyeDiagram{
	double period; 
	double freq;
	int window; // should be 2*period
	int sampleRate; 
	
	
	//The eye diagram consists of two palpebra (eyelids)
	//The supperior palbebra:
	//     _________
	//      \     /
	//       \___/
	vector< vector<node*> > palpebraSuperior;
	double* maxSuperior;	//the higher and lower envlope of the superior eyelid
	double* minSuperior;
	
	//The inferior palbebra:
	//         _____
	//        /     \
	//       /_______\
	//The inferior signal, models the transition from 1->1 / 1->1 / 1->0 / 0->1
	//Models the digital 
	vector< vector<node*> >  palpebraInferior;
	double* maxInferior;	// the higher and lower envlope of the inferior eyelid
	double* minInferior;

	//integration values
	//sum under signal one
	//sum over signal one

	//sum under signal two
	//sum over signal two

	//lebesgue 1
	//lebesgue 2
	//lebesgue 3
	//lebesgue 4

	public:
		EyeDiagram(double freq, int sampleRate, int w);
		~EyeDiagram();
		void push();
};