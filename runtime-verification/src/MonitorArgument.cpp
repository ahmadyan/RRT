#include "MonitorArgument.h"


MonitorArgument::MonitorArgument(int i){
	index = i ;
	monitorArgumentType = MONITOR_ARGUMENT_SIGNAL ;
}

MonitorArgument::MonitorArgument(double d){
	value = d ;
	monitorArgumentType = MONITOR_ARGUMENT_STATIC_TYPE ;
}

MonitorArgument::MonitorArgument(Monitor* m){
	monitor = m ;
	monitorArgumentType = MONITOR_ARGUMENT_DYNAMIC_TYPE ;
}

MonitorArgument::~MonitorArgument(void){}

int MonitorArgument::getIndex(){
	return index ;
}

double MonitorArgument::getValue(){
	if(monitorArgumentType == MONITOR_ARGUMENT_STATIC_TYPE ) return value ;
	//else if (monitorArgumentType == MONITOR_ARGUMENT_DYNAMIC_TYPE ) return monitor->getValue() ; //bad practise
	//else if (monitorArgumentType == MONITOR_ARGUMENT_SIGNAL ) return signal->getValue(index); 
}

Monitor* MonitorArgument::getMonitor(){
	return monitor ;
}

MonitorArgumentType MonitorArgument::getType(){
	return monitorArgumentType ;
}

