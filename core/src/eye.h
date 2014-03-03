#pragma once

#include <iostream>
#include <vector>

#include "config.h"
#include "node.h" 

using namespace std;
enum Transition { tnull, t01, t00, t10, t11, tboot };
class EyeDiagram{
	Configuration* config;

	double period; 
	double freq;
	double window; // should be 2*period
	double sampleRate; 
	int size;	//window/sampleRate

	int voltage;	//the index of the voltage variable that we are tracking
	
	//The eye diagram consists of two palpebra (eyelids)
	//The supperior palbebra:
	//     _________
	//      \     /
	//       \___/
	vector< vector<node*> > palpebraSuperior;
	double* maxSuperior;	//the higher and lower envlope of the superior eyelid
	double* minSuperior;
	int* maxSuperiorIndex; 
	int* minSuperiorIndex;

	//The inferior palbebra:
	//         _____
	//        /     \
	//       /_______\
	//The inferior signal, models the transition from 1->1 / 1->1 / 1->0 / 0->1
	//Models the digital 
	vector< vector<node*> >  palpebraInferior;
	double* maxInferior;	// the higher and lower envlope of the inferior eyelid
	double* minInferior;
	int* maxInferiorIndex;
	int* minInferiorIndex;

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
		EyeDiagram(Configuration* config);
		~EyeDiagram();
		void push(node*);
		void push(node* v, Transition tran);
		string toString();
		void test();
		void sum();

		int getWindowSize();

		node* getNode(int i);
};