#include "property.h"
#include "analogprop.h"

Property::Property(){}
Property::Property(string formulae){
	id = generateID() ;
	cout << "Property created " << id << endl ;
}


int Property::objectCount = 0;

int Property::generateID(){
	return objectCount++;
}

int Property::getID(){
	return id;
}

//creates an instance of the property populated with the node data
Property::Property(Property* prop, node* v){
	id = generateID();
	cout << "Property created " << id << endl;

	if (prop->argumentType == PTypeAnalog){

		cout << " ***** FOR DEBUGGING PURPOSE ONLY ***** " << endl;
		AnalogProperty* ap = (AnalogProperty*)prop; 
		Property* arg = prop->argument;
		AnalogProperty* ap2 = (AnalogProperty*)(prop->argument);

		argument = new AnalogProperty((AnalogProperty*)(prop->argument), v);
	}else if(argumentType == PTypeBoolean){
		//todo
	}
	argumentType = prop->argumentType;
}



double Property::eval(){
	if (argumentType == PTypeAnalog){
		return ((AnalogProperty*)argument)->eval();
	}
	else if (argumentType == PTypeBoolean){
		return argument->eval();
	}
	else{
		return argument->eval();
	}

}

/*
AnalogProperty(AnalogPropertyType, AnalogProperty* arg1);
AnalogProperty(AnalogPropertyType, AnalogProperty* arg1, AnalogProperty* arg2);
AnalogProperty(AnalogPropertyType, double);
AnalogProperty(AnalogPropertyType, int);
*/
Property::Property(int propID){
	id = generateID();
	if (propID == 1){	//property 1 for tdo circuit
		const int tdo_vd = 0;
		const int tdo_il = 1; 
		const int tdo_t = 2;

		AnalogProperty* ap_access_vc = new AnalogProperty(APT_ACCESS, tdo_vd);	//access vc
		AnalogProperty* ap_access_il = new AnalogProperty(APT_ACCESS, tdo_il);	//access iL
		AnalogProperty* ap_access_t = new AnalogProperty(APT_ACCESS, tdo_t);	//access time

		AnalogProperty* ap_const_03v = new AnalogProperty(APT_CONST, 0.3);		//constant 0.3 voltage
		AnalogProperty* ap_const_06v = new AnalogProperty(APT_CONST, 0.6);		//constant 0.6 voltage
		AnalogProperty* ap_const_003i = new AnalogProperty(APT_CONST, 0.03);		//constant 0.03 current
		AnalogProperty* ap_const_006i = new AnalogProperty(APT_CONST, 0.06);		//constant 0.06 voltage

		AnalogProperty* ap1 = new AnalogProperty(APT_LESS, ap_access_vc, ap_const_06v);
		AnalogProperty* ap2 = new AnalogProperty(APT_MORE, ap_access_vc, ap_const_03v);
		AnalogProperty* ap3 = new AnalogProperty(APT_LESS, ap_access_il, ap_const_006i);
		AnalogProperty* ap4 = new AnalogProperty(APT_MORE, ap_access_il, ap_const_003i);

		argument = ap1;
		argumentType = PTypeAnalog;
	}
	else if (propID == 2){
		const int tdo_vd = 0;
		const int tdo_il = 1;
		const int tdo_t = 2;

		AnalogProperty* ap_access_vc = new AnalogProperty(APT_ACCESS, tdo_vd);	//access vc
		AnalogProperty* ap_access_il = new AnalogProperty(APT_ACCESS, tdo_il);	//access iL
		AnalogProperty* ap_access_t = new AnalogProperty(APT_ACCESS, tdo_t);	//access time

		AnalogProperty* ap_const_03v = new AnalogProperty(APT_CONST, 0.3);		//constant 0.3 voltage
		AnalogProperty* ap_const_06v = new AnalogProperty(APT_CONST, 0.6);		//constant 0.6 voltage
		AnalogProperty* ap_const_003i = new AnalogProperty(APT_CONST, 0.03);		//constant 0.03 current
		AnalogProperty* ap_const_006i = new AnalogProperty(APT_CONST, 0.06);		//constant 0.06 voltage

		AnalogProperty* ap2 = new AnalogProperty(APT_MORE, ap_access_vc, ap_const_03v);

		argument = ap2;
		argumentType = PTypeAnalog;
	}
	else if (propID == 3){
		const int tdo_vd = 0;
		const int tdo_il = 1;
		const int tdo_t = 2;

		AnalogProperty* ap_access_vc = new AnalogProperty(APT_ACCESS, tdo_vd);	//access vc
		AnalogProperty* ap_access_il = new AnalogProperty(APT_ACCESS, tdo_il);	//access iL
		AnalogProperty* ap_access_t = new AnalogProperty(APT_ACCESS, tdo_t);	//access time

		AnalogProperty* ap_const_03v = new AnalogProperty(APT_CONST, 0.3);		//constant 0.3 voltage
		AnalogProperty* ap_const_06v = new AnalogProperty(APT_CONST, 0.6);		//constant 0.6 voltage
		AnalogProperty* ap_const_003i = new AnalogProperty(APT_CONST, 0.03);		//constant 0.03 current
		AnalogProperty* ap_const_006i = new AnalogProperty(APT_CONST, 0.06);		//constant 0.06 voltage

		AnalogProperty* ap3 = new AnalogProperty(APT_LESS, ap_access_il, ap_const_006i);
		
		argument = ap3;
		argumentType = PTypeAnalog;
	}
	else if (propID == 4){
		const int tdo_vd = 0;
		const int tdo_il = 1;
		const int tdo_t = 2;

		AnalogProperty* ap_access_vc = new AnalogProperty(APT_ACCESS, tdo_vd);	//access vc
		AnalogProperty* ap_access_il = new AnalogProperty(APT_ACCESS, tdo_il);	//access iL
		AnalogProperty* ap_access_t = new AnalogProperty(APT_ACCESS, tdo_t);	//access time

		AnalogProperty* ap_const_03v = new AnalogProperty(APT_CONST, 0.3);		//constant 0.3 voltage
		AnalogProperty* ap_const_06v = new AnalogProperty(APT_CONST, 0.6);		//constant 0.6 voltage
		AnalogProperty* ap_const_003i = new AnalogProperty(APT_CONST, 0.03);		//constant 0.03 current
		AnalogProperty* ap_const_006i = new AnalogProperty(APT_CONST, 0.06);		//constant 0.06 voltage

		AnalogProperty* ap4 = new AnalogProperty(APT_MORE, ap_access_il, ap_const_003i);
		argument = ap4;
		argumentType = PTypeAnalog;
	}
	cout << "Property created " << id << endl ;
}
