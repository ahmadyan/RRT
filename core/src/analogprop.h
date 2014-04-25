#pragma once

#include <iostream>
#include <string>
#include "property.h"
#include "node.h"
#include "deviation.h"

using namespace std;

enum AnalogPropertyType {APT_ACCESS, APT_PLUS, APT_MINUS, APT_MULT, APT_DIV, APT_NORM, APT_DIST, APT_DEV, APT_EQ, APT_LESS, APT_MORE, APT_CONST };
class AnalogProperty : public Property{	
	node* v;	//the data itself
	int var;	//which dimension in v we are reasoning with? 
	bool executed;	//used in eval, we execute each property only once for the incremental operation
	const double epsilon = 1e-6;

	//for efficient deviation implementation, we create a seperate class to hold the signal evnelope called deviation
	Deviation* deviation;
	//because there are arbitrary number of signal envelopes, the deviationID refers to each of them uniquely. 
	int deviationID;

	AnalogPropertyType op;
	vector<AnalogProperty*> arguments; 
	double result;

	public:
	//Construct the property from given formulae, typically used only once at the begining of the algorithm
	AnalogProperty(AnalogPropertyType, AnalogProperty* arg1);
	AnalogProperty(AnalogPropertyType, AnalogProperty* arg1, AnalogProperty* arg2);
	AnalogProperty(AnalogPropertyType, double);
	AnalogProperty(AnalogPropertyType, int);

	//Copy constructor, construct the property from another instance, without evaluating it. 
	//Used at every iteration to put a new property parse tree at every node
	AnalogProperty(AnalogProperty*, node*);
	AnalogProperty();
	~AnalogProperty();

	void registerDeviationClass(Deviation*);
	void setDeviationID(int);
	void setVariable(int);

	double eval();

	void setNode(node*); 
	node* getNode();
};