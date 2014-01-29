#include "RRT.h"
#include <iostream>
#include <limits>
#include <stack>
using namespace std;

RRT::RRT(int _d, int _k, int _var, string nam){
    d=_d;
    k=_k;
	var=_var;
    name=nam;
    min = new double[d];
    max = new double[d];
	variationMin = new double[var];
	variationMax = new double[var];
    //default values for min-max
    for(int i=0;i<d;i++){
        min[i]=0;
        max[i]=1;
    }
	for (int i = 0; i < var; i++){
		variationMin[i] = 0;
		variationMax[i] = 0;
	}
}

RRT::RRT(string fileName){
	load(fileName);
}

RRT::~RRT(){
    //if(min!=NULL) delete min ;
    //if(max!=NULL) delete max ;
}


double RRT::unifRand(){
    return rand() / double(RAND_MAX);
}
double RRT::unifRand(double a, double b){
    return (b-a)*unifRand() + a;
}

node* RRT::getNearestNode(node* q_sample){
	return root->getNearestNode(q_sample, min, max, false).first ;
}

//There is no actual simulation in this code,
//except a pretty cool operator overloading!
void RRT::buildUniform(double* initialState){
    root = new node(d);
    root->set(initialState);
	root->setRoot();
    for(int i=0;i<k; i++){
        //create a new sample
        node* q_sample = new node(d);
        q_sample->randomize(min, max);
		node* q_near = getNearestNode(q_sample);
        //find the best trajectory using integration or simulation, constructing the new node q_new
        double delta = 0.5;
        node q_new = *q_near + ( ( *q_sample - *q_near ) * delta ) ;
        q_new.dump();
		q_near->addChildren( &q_new ) ;
		q_new.setParent(q_near);
    }
}

void RRT::build(double* initialState){
    root = new node(d);
    root->set(initialState);
	root->setRoot();
	
    for(int i=0;i<k; i++){
        node* q_sample = new node(d);
        q_sample->randomize(min, max);
        node* q_near = getNearestNode(q_sample);
    
		double* state_near = q_near->get() ;
        double* ic = new double[d];
        for(int j=0;j<d;j++) ic[j]=state_near[j];
        double delta = 0.5;
        //variation or input to the system
        //todo: should be defined in the main or system, not here
		vector<double> param ;
		for (int i = 0; i < var; i++){
			param.push_back(unifRand(0.29, 0.31));
		}
		vector<string> settings;
		vector<double> result = system->simulate(ic, param, settings, delta);

		double* state = new double[d];
		for (int i = 0; i < d; i++){
			state[i] = result[i+1];
		}
		delete ic;
        cout << "state*=" << state[0] << " " << state[1] << endl ;
        node* q_new = new node(d);
        q_new->set(state);

        q_near->addChildren( q_new ) ;
		q_new->setParent(q_near);
	}
}


//Define the minimum and maximum value for each dimension
void RRT::setBound(int i, double _min, double _max){
    min[i] = _min;
    max[i] = _max;
}

void RRT::setVariationBound(int i, double _min, double _max){
	variationMin[i] = _min;
	variationMax[i] = _max;
}

/*
node RRT::getNode(int i){
    return nodes[i].second;
}

int RRT::getNodeParentId(int i){
	return nodes[ nodes[i].first ].first;
}

//returns the parent of node i in the tree, we hold the parent# in nodes[i].first
node RRT::getParentNode(int i){
    cout << "get parent for node " << i << " = " << nodes[i].first << endl ;
    return nodes[ nodes[i].first ].second;
}

int RRT::sampleNumber(){
    return nodes.size();
}
*/
int RRT::getDimension(){
    return d; 
}

void RRT::setSystem(System* s){
    system = s ;
}

double RRT::getMin(int i){
    return min[i] ;
}

double RRT::getMax(int i){
    return max[i];
}

//Saving the RRT into a file
//Stores each-node-id, parent-node-id, node-data, input
void RRT::save(string fileName){
	cout << "Saving the RRT" << endl ;
	ofstream file;
	file.open (fileName.c_str());
	file << "rrt" << endl;
	file << d << endl ;		//number of dimensions
	file << var << endl ;
	file << k << endl ;		//number of nodes

	//saving the bounds on each dimensions
	for(int i=0;i<d; i++){
		file << min[i] << endl ;
		file << max[i] << endl ;
	}

	//start from the root and recursively print every node
	//root->save(file);

	//linear printing of every nodes
	//nodes are being printed in this format: id parentid data timestamp

	cout << "Nodes.size="<< nodes.size() << endl ;
	for(int i=0;i<nodes.size();i++){
		cout << "Saving node " << i << endl ;
		file << i << " ";

		if(nodes[i]->isRoot()){
			file << "-1" << " " ;
		}else{
		for(int j=0;j<nodes.size();j++){
			if(	nodes[i]->getParent()->getID() == nodes[j]->getID() )
				file << j << " " ;
			}
		}
		
		for(int j=0;j<nodes[i]->getDimension();j++){
			file << nodes[i]->get(j) << " " ;
		}
		file << endl ;
	}

	file.close();
}

void RRT::load(string fileName){
	string line;
	ifstream file (fileName.c_str());
	if (file.is_open()){
		file >> name ;
		cout << "Loading RRT " << name << endl ;
		file >> d ;
		file >> var ;
		file >> k;

		min = new double[d];
		max = new double[d];
		cout << "d,k="<< d  << " " << k << endl ;
		for(int i=0;i<d; i++){
			file >> min[i] ;
			file >> max[i] ;
			cout << min[i] << " " << max[i] << endl ;
		}
		cout << "loading nodes" << endl ;
		

		node* lastParent; 

		for(int i=0;i<k;i++){
			double* data = new double[d];
			int id =-1;
			int parent_id = -1 ;
			file >> id ;
			file >> parent_id ;
			for(int j=0;j<d;j++){
				file >> data[j];
			}

			cout << "Loading node " << i << endl ;
			node newNode= node(d, id, data);
			nodes.push_back(&newNode);
			

			//just adding the parent is not enought, we need a proper get node by id method. 





			cout << newNode.toString() << endl ;
		}
		file.close();
	}
}

string RRT::getName(){
	return name ;
}

void RRT::setdt(double d){
	dt = d ;
}

double RRT::getdt(){
	return dt;
}

node* RRT::getRoot(){
	return root;
}


vector<node*> RRT::getNearestNode(node* q_sample, double errorTolerance, bool timed){
	//old method, O(n), used recursive tree search
	//return root->getNearestNode(q_sample, min, max, true).first ;

	vector<node*> results;
	if(errorTolerance <= 0){ 
		//Searching for the nearest node, 
		//Usually used for searching for closest node in the RRT from q_sample during the RRT loop
		double closestDistance= ( timed ? q_sample->timed_distance(nodes[0], min, max) : q_sample->distance(nodes[0], min, max) );
		int closestNodeIndex=0; 
		for(int i=0; i<nodes.size();i++){
			double d= ( timed ? q_sample->timed_distance(nodes[i], min, max) : q_sample->distance(nodes[i], min, max) );
			if( d<closestDistance ){
				closestDistance=d;
				closestNodeIndex=i;
			}
		}
		results.push_back(nodes[closestNodeIndex]);
	}else{
		//There is a positive errorTolerance, which means we are searching for the set of close nodes within the errorTolerance epsilon
		//THis is usually used for the casting searches (i.e. nodes that are within the same state, but possibly different time)
		for(int i=0; i<nodes.size();i++){
			double d= ( timed ? q_sample->timed_distance(nodes[i], min, max) : q_sample->distance(nodes[i], min, max) );
			cout << "Searching for closest neighbors" << endl ;
			cout << "distance=" << d << endl ;
			if(d<=errorTolerance){
				if(nodes[i]->getID()!=q_sample->getID()){
					results.push_back(nodes[i]);
				}
			}
		}
	}
	return results;
}


void RRT::addMonitor(Monitor* m){
	cout << "Adding a new monitor to the RRT , total monitors=" << monitors.size() << endl ;
	monitors.push_back(m);
}
