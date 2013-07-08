#pragma once

using namespace std;

class Property{
	
	int propertyID;				//unique
    static int propertyCount ;
	
	//methods
	static int generateMonitorID();
public:
    Property();
    ~Property();
    bool isSatisfied();
	int getPropertyID();
};
