#pragma once

#include <vector>
#include <utility>
#include "node.h"
#include "System.h"

using namespace std;

class RRT{
protected:
    int k ; // total number of nodes 
    int d ; // tree dimension, each node has a d dimension. (d-1 state + 1 time)
    //vector< pair<int,node> > nodes; // this is the simple tree structure, each pair consists of the id to the parent and the node itself
	node* root ;
    double* min ;
    double* max ;
    System* system ;
    string name;
	double dt ;
public:
    RRT(int, int, string);
    RRT(string);
    ~RRT();
    void load(string fileName);
    void save(string fileName);
    string getName();
    double unifRand();
	double unifRand(double a, double b);
    void build(double*);
    void buildUniform(double* initialState);
    void setBound(int d, double _min, double _max);
    //int sampleNumber();
    //node getNode(int i);
    //node getParentNode(int i);
	//int getNodeParentId(int i);
    int getDimension();
    double getMin(int i);
    double getMax(int i);
    void setSystem(System*);
	void setdt(double);
	double RRT::getdt();
	node* getNearestNode(node* q_sample);
	node* getRoot();
};
