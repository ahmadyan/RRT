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
    //default values for min-max
    for(int i=0;i<d;i++){
        min[i]=0;
        max[i]=1;
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
        double* state = new double[d];
        for(int j=0;j<d;j++) state[j]=state_near[j];
        double delta = 0.5;
        //variation or input to the system
        //todo: should be defined in the main or system, not here
		double* param = new double[var] ;
        param[0] = unifRand(0.29, 0.31);
        //double param = 1.4 ;
        cout << "state==" << state[0] << " " << state[1] << endl ;
        system->simulate(state, param, delta);
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
	root->save(file);
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
		int dummy ;
		for(int i=0;i<k;i++){
			double* data = new double[d];
			int parent_id = 0 ;
			file >> dummy ;
			file >> parent_id ;
			for(int j=0;j<d;j++){
				file >> data[j];
			}
			node newNode = node(d);
			newNode.set(data);

			

			nodes.push_back( make_pair(parent_id, &newNode));
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