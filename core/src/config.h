//This is a code from Daedalus.

#pragma once
#include <map>
#include <string>
using namespace std;

//#pragma warning( disable : 4786)
//#pragma warning( disable : 4018)
//#pragma warning( disable : 4224)
//#pragma warning( disable : 4244)
//#pragma warning( disable : 4996)
//#pragma warning( disable : 4291) //no matching operator delete
//#pragma warning( disable : 4305) // truncation from 'const double' to 'float'

class Configuration{
    map<string, string> db;

public:
    Configuration(string);
	string get(string);
    void getParameter(string, string*);
    void getParameter(string, double*);
	void getParameter(string parameter, int i, string* result);
	void getParameter(string, int, double*);
    void getParameter(string parameter, int* result);
    void setParameter(string, string);
	bool checkParameter(string key, string value);
	bool checkParameter(string parameter, int i, string value);
};
