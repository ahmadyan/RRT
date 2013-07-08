#pragma once
#include <vector>
#include <map>
#include "System.h"
#include "Rrt.h"

//#include "MonitorArgument.h"
using namespace std;

class Monitor{
	Property* property; //Every monitor is checking a unique property associated with it during construction

	int monitorID;				//unique
    static int monitorCount ;
	
	//methods
	static int generateMonitorID();

public:
	Monitor(Property* p);
	~Monitor();
    int getMonitorID();

	int monitor_update_add_state();
	int monitor_update_add_edge();
	int monitor_get_solution();
};



//#define MONITOR_EXCEPTION_DIVISION_BY_ZERO	0x142495E
//#define MONITOR_EXCEPTION_VIOLATION			0x0F588C9
//#define MONITOR_EXCEPTION					0x0B8F500
//#define epsilon								1e-9
//	This enum type will contain the different type of monitor that our algorithm supports, 
//	both analog and temporal. In our implementation, we don't differentiate between analog and temporal.
/*enum MonitorType { 
	CONSTANT, 
	ANALOG_BINARY_ADD,
	ANALOG_BINARY_SUB,
	ANALOG_BINARY_MUL,
	ANALOG_BINARY_DIV,
	ANALOG_SHIFT,
	ANALOG_NORM_L2,
	
	ANALOG_MIN_SIBLING,
	ANALOG_MAX_SIBLING,
	ANALOG_DIFF_SIBLING,
	//logical operators:
	LOGIC_LESS_THAN, 
	LOGIC_GREATER_THAN,
	LOGIC_EQUAL,
	LOGIC_LESS_THAN_OR_EQUAL,
	LOGIC_GREATER_THAN_OR_EQUAL,
	LOGIC_AND,
	LOGIC_OR,
	LOGIC_EVENTUALLY,
	LOGIC_UNTIL,
	LOGIC_JITTER
};

enum MonitorStatus {DONTKNOW, SATISFIED, VIOLATED};
*/
/*
class Monitor{
	RRT* rrt ;
	
	vector<MonitorArgument*> arguments ;
	MonitorType type ;

	//variables used for model checking
	map<node*, int> track ;	//Holds number of Phi till now
	map<node*, bool> trend ; //Holds if Phi was always true till now
	map<node*, MonitorStatus> status ;
	map<node*, bool> updated ;
public:
    Monitor(RRT* rrt, MonitorType type);
	void push_back(double data);
	void Monitor::push_back(int data);
	void push_back(Monitor*);
	~Monitor();
	static int generateMonitorID();
    int getMonitorID();
	//This is the main function, this function will be called after adding each new node to the RRT
	double verify(node* q_new);
	double getValue();
	double getArgumentValue(int i, node* q);
	bool analyze(int i, node* q_new);
	void update(node* q_new, bool value);
	void update(node* q);
};*/