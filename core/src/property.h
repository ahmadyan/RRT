#pragma once

#include <iostream>
#include <string>
#include "node.h"

using namespace std;

enum PropertyArgumentType	{PTypeAbstract, PTypeAnalog, PTypeBoolean};
class Property{
public:
	static int objectCount ;
	int id;
	int result;

	Property* argument; 
	


	PropertyArgumentType argumentType;
	//Construct the property from given formulae, typically used only once at the begining of the algorithm
	Property();
	Property(string);		
	Property(int);	//same as above, but hard-coded
	Property(Property*, node*);
	
	static int generateID();
	int getID();

	double eval();
};