#pragma once

#include <iostream>
#include <string>
#include "node.h"

using namespace std;

class Property{
	static int objectCount ;
	int id;
	node* v;
	//vector<ops>

	int result;

	public:
	//Construct the property from given formulae, typically used only once at the begining of the algorithm
	Property(string);		
	Property(int);	//same as above, but hard-coded

	//Copy constructor, construct the property from another instance, without evaluating it. 
	//Used at every iteration to put a new property parse tree at every node
	Property(Property*);	
	~Property();
	static int generateID();
	int getID();

	void setNode(node*);
	node* getNode();

	int getResult();

	void check();
};