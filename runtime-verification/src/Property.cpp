#include "Property.h"
int Property::propertyCount = 0 ;
Property::Property(){
	propertyID = generateMonitorID() ;
}

Property::~Property(){
}

bool Property::isSatisfied(){
    return true ;
}

int Property::generateMonitorID(){
	return propertyCount++;
}  

int Property::getPropertyID(){
	return propertyID;    
}