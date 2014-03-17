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
	
	int vSize;	//for lebesge integrals dv
	double dv;

	int voltage;	//the index of the voltage variable that we are tracking

	double* nadir;

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
	double* leftSuperior;
	double* rightSuperior;
	int* leftSuperiorIndex;
	int* rightSuperiorIndex;

		

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
	double* leftInferior;
	double* rightInferior;
	int* leftInferiorIndex;
	int* rightInferiorIndex;


	double maxOne;
	double minOne;
	double maxZero;
	double minZero;

	//these two sets contains the nodes that we can initiate a jitter from, usually nodes with t<t_jitter_max
	vector<node*> jitterFrontierSet10; 
	vector<node*> jitterFrontierSet01;

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