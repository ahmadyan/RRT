#pragma once
#include <iostream>
#include <vector>
using namespace std;

class transient{
    int n ;//dimension
    vector<double> start ;  //x_0, x_1, ..., x_n
    vector<double> end   ;  //y_0, y_1, ..., y_n
    double t ;
    double param ;
    
public:
    transient(int _n, double _t, double _param);
    int getN();
    double getT();
    double getParam();
    vector<double> getStart();
    vector<double> getEnd();
    double getStart(int i);
    double getEnd(int i);
    void dump();
    void setStartingPoint(double* point);
    void setEndPoint(double* point);

};