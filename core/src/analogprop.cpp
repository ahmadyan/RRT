#include "analogprop.h"
AnalogProperty::AnalogProperty(){}

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
	result = d;
}

AnalogProperty::AnalogProperty(AnalogPropertyType t, int d){
	op = t;
	var = d;
}

AnalogProperty::AnalogProperty(AnalogProperty* property, node* nd){
	cout << "Copy constructing node for " << nd->getID() << endl;
	op = property->op;
	var = property->var;
	deviationID = property->deviationID;
	cout << " argumentSize = " << property->arguments.size() << endl;
	v = nd;
	for (int i = 0; i < property->arguments.size(); i++){
		//this line will also copy construct every property in the arguments recursively.
		//so in order for this to work, we have to copy only the root analog property, not every argument seperately. 
		arguments.push_back(new AnalogProperty(property->arguments[i], nd));
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
			for (int i = 0; i < v->getDimension()-1; i++){
				result += v->get(i)*v->get(i);
			}
			result = sqrt(result);
			break;  
		case APT_DIST: 
			//todo: implement this. 
			for (int i = 0; i < v->getDimension() - 1; i++){
				result += v->get(i)*v->getParent()->get(i);
			}
			result = sqrt(result);
			break;  
		case APT_DEV:
			//executed = false;
			if (var == -1){
				//compute the deviation of the result
				result = arguments[0]->eval();
			}
			else{
				result = v->get(var);
			}
			deviation->push(v->getTime(), result, 0);
			result = deviation->get(v->getTime(), 0);
			 
			break;  
	
		//--Comparison Operators---------------------------------
		case APT_EQ: 
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
	cout << "Evaluation of Property " << op << " is " << result << endl;
	return result;
}

void AnalogProperty::registerDeviationClass(Deviation* d){
	deviation = d;
}

void AnalogProperty::setDeviationID(int x){
	deviationID = x;
}

void AnalogProperty::setVariable(int x){
	var = x;
}