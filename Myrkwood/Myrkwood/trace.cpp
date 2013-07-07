#include "trace.h"

Trace::Trace(){
}
Trace::~Trace(){
}

void Trace::AddSample(vector<double> sample){
    data.push_back(sample);
}

vector<double> Trace::getSample(int i){
    return data.at(i);
}

int Trace::getSampleSize(){
    return data.size();
}

vector<vector<double> > Trace::getData(){
    return data;
}