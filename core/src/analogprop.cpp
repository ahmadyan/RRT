#include "analogprop.h"


AnalogProperty::AnalogProperty(AnalogPropertyType t, AnalogProperty* arg1){
	op = t;
	arguments.push_back(arg1);
	executed = false;
}
AnalogProperty::AnalogProperty(AnalogPropertyType t, AnalogProperty* arg1, AnalogProperty* arg2){
	op = t;
	arguments.push_back(arg1);
	arguments.push_back(arg2);
	executed = false;
}
AnalogProperty::AnalogProperty(AnalogPropertyType t, double d){
	op = t;
	if (t == APT_CONST){
		result = d;
		executed = true;
	}
}

AnalogProperty::AnalogProperty(AnalogPropertyType t, int d){
	op = t;
	if (t == APT_ACCESS){
		var = d;
		executed = true;
	}
}

AnalogProperty::AnalogProperty(AnalogProperty* property){
	op = property->op;
	var = property->var;
	for (int i = 0; i < property->arguments.size(); i++){
		//this line will also copy construct every property in the arguments recursively.
		//so in order for this to work, we have to copy only the root analog property, not every argument seperately. 
		arguments.push_back(new AnalogProperty(property->arguments[i]));
	}
	executed = false;
	//we do not copy the node itself

	if (op == APT_CONST){
		executed = true;
		result = property->result;
	}
}

AnalogProperty::~AnalogProperty(){
}

void AnalogProperty::setNode(node* _v){
	v=_v;
}

node* AnalogProperty::getNode(){
	return v;
}

double AnalogProperty::eval(){
	if (executed) return result;
	result = 0;
	executed = true;
	switch (op){
		//-----------------------------------------------------
		case APT_ACCESS:
			result = v->get(var);
			break;

		case APT_CONST:
			//do nothing, already addressed in constructor
			break;

		//--Mathematical Operators------------------------------
		case APT_PLUS:
			result = arguments[0]->eval() + arguments[1]->eval();
			break;
		case APT_MINUS: 
			result = arguments[0]->eval() - arguments[1]->eval();
			break; 
		case  APT_MULT: 
			result = arguments[0]->eval() * arguments[1]->eval();
			break;  
		case APT_DIV: 
			result = arguments[0]->eval() / arguments[1]->eval();
			break;  

		//--Distance Operators---------------------------------
		case APT_NORM: 
			result = arguments[0]->eval() + arguments[1]->eval();
			break;  
		case APT_DIST: 
			result = arguments[0]->eval() + arguments[1]->eval();
			break;  
		case APT_DEV: 
			result = arguments[0]->eval() + arguments[1]->eval();
			break;  
	
		
		//--Comparison Operators---------------------------------
		case APT_EQ: 
			double epsilon = 1e-6;
			if (abs( arguments[0]->eval() - arguments[1]->eval() ) < epsilon){
				result = 1;
			}
			else{
				result = 0;
			}
			break;  
		case APT_LESS: 
			if (arguments[0]->eval() < arguments[1]->eval()){
				result = 1;
			}
			else{
				result = 0;
			}
			break;  
		case APT_MORE: 
			if (arguments[0]->eval() > arguments[1]->eval()){
				result = 1;
			}
			else{
				result = 0;
			}break;

		default:
			cout << "Uknown analog operators: " << op << endl;
			exit(1);
			break;
	}
	return result;
}