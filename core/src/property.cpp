#include "property.h"
#include "pll.h"

Property::Property(string formulae){
	id = generateID() ;
	cout << "Property created " << id << endl ;
}

Property::Property(int){
	id = generateID() ;
	cout << "Property created " << id << endl ;
}

Property::Property(Property* property){
	id = generateID() ;
	cout << "Property constructed " << id << endl ;
}

Property::~Property(){
}

int Property::objectCount = 0 ;

int Property::generateID(){
	return objectCount++;
}

int Property::getID(){
	return id;
}

void Property::setNode(node* _v){
	v=_v;
}

node* Property::getNode(){
	return v;
}

int Property::getResult(){
	return result;
}

void Property::check(){
	//debug: a very simple check
	cout << endl << endl << endl << "---------------------------------" << endl ;
	double t=v->get(pll_time);
	double diffVoltage = abs( v->get(pll_e) - v->get(pll_eb)) ;
	if( t>1e-5 ){
		if( diffVoltage>0.03 ){
			result = 0;
		}else{
			result = 1;
		}
	}

	cout << v->get(pll_time) << endl ;
	
	cout << diffVoltage << endl;
	cout << "RESULT=" << result << endl ;
	cout << endl << endl << endl << "---------------------------------" << endl ;
	
}



/*
	Checking safety
	Phi = \nu z. (q ^ next z)
*/