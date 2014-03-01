#include "transient.h"

transient::transient(int _n, double _t, double _param){
    n=_n;
    t=_t;
    param=_param ;
}

void transient::setStartingPoint(double* point){
    for(int i=0;i<n;i++)
        start.push_back(point[i]);
}

void transient::setEndPoint(double* point){
    for(int i=0;i<n;i++)
        end.push_back(point[i]);
}

int transient::getN(){
    return n;
}

double transient::getT(){
    return t;
}

double transient::getParam(){
    return param ;
}

vector<double> transient::getStart(){
    return start;
}

vector<double> transient::getEnd(){
    return end;
}

double transient::getStart(int i){
    return start[i];
}

double transient::getEnd(int i){
    return end[i];
}
void transient::dump(){
    cout << "This transient is from " ;
    for(int i=0;i<n;i++) cout << start[i] << " " ;
    cout << " -->to--> " ;
    for(int i=0;i<n;i++) cout << end[i] << " " ;
    cout << " for " << t <<  " and param=" << param << endl ;
}