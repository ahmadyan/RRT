#pragma once

#include <vector>
#include <utility>
#include "node.h"
#include "System.h"
#include "Monitor.h"
#include "config.h"
#include "eye.h"
#include "config.h"
using namespace std;
class Monitor;

class RRT{
protected:
	Configuration* config;
    int k ; // total number of nodes 
    int d ; // tree dimension, each node has a d dimension. (d-1 state + 1 time)
    vector<node*> nodes; //keeps a copy of each node for direct access, pair<parent_id, child>
	node* root ;
    double* min ;
    double* max ;
	double* variationMin;
	double* variationMax;
    System* system ;
    string name;
	double dt ;

	int var; //number of variation parameters
	EyeDiagram* eye;

	vector<Monitor*> monitors;
public:
    RRT(Configuration*, int, int, int, string);
	RRT(Configuration* c, string);
    ~RRT();
    void load(string fileName);
    void save(string fileName);
    string getName();
    
	double generateNormalSample(double mean, double std);
	double generateTruncatedNormalSample(double mean, double std, double min, double max);
	double generateUniformSample();
	double generateUniformSample(double a, double b);

    void build(double*);
    void buildUniform(double* initialState);
    void setBound(int d, double _min, double _max);
	void setVariationBound(int i, double _min, double _max);
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


	vector<node*> getLatestNode();
	vector<node*> NearestNodeInProjectiveSpace(node* q_sample);
	vector<node*> NearestNodeUsingTimeDistance(node* q_sample);
	vector<node*> NearestNodeUsingEuclideanDistance(node* q_sample);
	vector<node*> recursiveNearestNodeSearch(node* q_sample);
	vector<node*> randomizedNearestNodeSearchEuclideanDistance(node* q_sample);
	vector<node*> randomizedNearestNodeSearchTimeDistance(node* q_sample);
	vector<node*> getNearestNode(node* q_sample);
	node* getRoot();
	node* getNode(int);
	void addMonitor(Monitor*);
	void setConfig(Configuration*); 

	EyeDiagram* getEyeDiagram();
};
