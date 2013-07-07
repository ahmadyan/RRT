#include "node.h"
#include <cmath>
#include <iostream>
using namespace std;

Node::Node(geometry::Point* p, bool r){
    point = new geometry::Point(p);
    root = r;
}

Node::~Node(){
    //delete point;
}

bool Node::isRoot() const{
    return root;
}

//compute the norm between current node and given node a
double Node::distance(Node* a, Circuit* circuit){
    return distance(a->point, circuit);
}

//Normalized distance, requires the circuit class.
double Node::distance(geometry::Point* p, Circuit* circuit){
	double d = 0;
	for(int i=0;i<circuit->getDimension();i++){
        double di = point->getData(i) - p->getData(i) ;
        double di_normalized = di / (circuit->getMax(i)-circuit->getMin(i));
		d +=  di_normalized*di_normalized;
	}
	d = sqrt(d);
	return d;
}

double* Node::getState(){
    return point->getData();
}

Node* Node::getParent(){
    return parent;
}

void Node::setParent(Node* p){
    parent = p ;
}

Type Node::getTreeType(){
    return treeType;
}

void Node::setTreeType(Type t){
    treeType=t;
}

Type Node::getNodeType(){
    return nodeType;
}

void Node::setNodeType(Type t){
    nodeType=t;
}

string Node::toString(){
    return point->toString();
}

geometry::Point* Node::getPoint() {
    return point ;
}

double Node::getData(int i) {
    return point->getData(i) ;
}

void Node::setTime(double t0){
    t=t0;
}

double Node::getTime(){
    return t;
}
