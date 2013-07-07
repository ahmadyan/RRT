#include "sampling.h"

RPG::RPG(Circuit* c){
    circuit = c;
}

RPG::~RPG() {
}


double RPG::random(double a, double b){
    return (b-a)*(rand() / double(RAND_MAX)) + a;
}


geometry::Point* RPG::sample(){
	geometry::Point* p = new geometry::Point(circuit->getDimension()) ;
    for(int i=0;i<circuit->getDimension(); i++){
        double x = random( circuit->getMin(i), circuit->getMax(i) );
        p->setData(i, x);
    }
	return p;
}
