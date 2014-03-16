#include "node.h"
#include <iostream>
#include <math.h>
#include <fstream>
#include <cstdlib>
using namespace std;

node::node(int _n){
	jitter = 0;
	id = generateID();
	root=false; 
	n=_n;
	data = new double[n];
	for(int i=0;i<n;i++) data[i]=0;
}

node::node(int _n, int _id, double* _data){
	id=_id; 
	root=false;
	n=_n;
	data=_data;
	if(id>=objectCount)
		objectCount=id+1;
}

//copy constructor
node::node(const node& Node){
	id = generateID() ;
	root = Node.root ;
	n = Node.n ;
	data = new double[n];
	for(int i=0;i<n;i++){
		data[i] = Node.data[i];
	}
	jitter = 0;
}

//This piece of code causes unknown problems with visual studio 2012.
node::~node(){
	//for(int i=0;i<children.size();i++){
	// delete children[i];
	//}
	//	delete data ;
}


int node::objectCount = 0 ;

int node::generateID(){
	return objectCount++;
}

int node::getID(){
	return id;
}

double node::unifRand(){
	return rand() / double(RAND_MAX);
}
double node::unifRand(double a, double b){
	return (b-a)*unifRand() + a;
}

#include <random>
double node::normalRand(double min, double max){
	std::random_device rd;
	std::mt19937 generator(rd());
	double mean = (min+max)/2;
	double std = 0.4;
	std::normal_distribution<double> normal(mean, std);

	double d=normal(generator);
	if (d < min) d = min;
	if (d>max) d = max;
	return d;
}

void node::set(double* _x){
	//todo: check for size x
	data=_x;
}

void node::set(int i, double d){
	data[i]=d;
}


/// This setter should only be called from the TimedRRT, not the RRT class
void node::set(vector<double> v){
	//First, copying the data
	for (int i = 0; i < n-1; i++){
		data[i] = v[i + 1];
	}
	//Then copying the time-stamp
	data[n - 1] = v[0];
}

double* node::get(){
	return data ;
}
double  node::get(int i){
	return data[i];
}

int node::getDimension(){
	return n;
}

double node::getTime(){
	return data[n-1];
}

//assigns a uniform random value to the new node,
//this methods is usually called for generatign a new sample in the RRT, so
//the RRT can find the nearest node in the tree to this newly generated node.
//random number is between (min[i], max[i]) range.
void node::randomize(double* min, double* max){
	for(int i=0;i<n;i++){
		data[i] = normalRand(min[i], max[i]);
		//data[i]=unifRand(min[i], max[i]);
	}
}

//compute the norm between current node and given node a
double node::distance(node* a, double* max, double* min){
	double d = 0;
	for(int i=0;i<n;i++){
		d +=  ( ( data[i]-a->data[i] ) * ( data[i] - a->data[i] ) / (max[i]-min[i]));
	}
	d = sqrt(d);
	return d;
}

//compute the norm between current node and given node a
//In timed_distance, unlike normal distance function, we add an extra emphasis on the time dimension.
//We wish to ensure forwawrd progress in time, therefore we are going to pick the samples that 
//are more forward in time (i.e. their time dimension is closer to the time_envlope, the latest time sampled discovered so far, 
//with higher probability. 
//How? We scale samples by 1-time_difference.
double node::timed_distance(node* a, double* max, double* min){
	double d = 0;
	double weightFactor = unifRand(0, 10);			//how much emphasize should we put it time dimension?
	double sampleTime = data[n - 1]; //time_envlope is the latest time discovered so far. 
	for (int i = 0; i < n - 1; i++){
		d += ((data[i] - a->data[i]) * (data[i] - a->data[i]) / ((max[i] - min[i])*(max[i] - min[i])));
	}
	if (sampleTime >= 0){
		double relativeTimeDistance = abs(a->data[n - 1] - sampleTime) / sampleTime;  // 0 <= pos <= 1, 0 is close in time where as 1 is faraway in time
		d += weightFactor*(relativeTimeDistance)*d;
		//d += weightFactor*(1 - relativeTimeDistance)*d;
	}
	return d;
}

//TODO: issue a warning when the vector size does not match in any of overloaded operators
node& node::operator=(const node &rhs){
	if(this==&rhs) return *this; 
	delete data ;   //deallocate previous copy of data, (we can's reuse it, since the n may be different between two node)
	n = rhs.n;
	data = new double[rhs.n];
	for(int i=0;i<n;i++){
		data[i] = rhs.data[i];
	}
	return *this;
}

node& node::operator+=(const node &rhs){
	for(int i=0;i<n;i++){
		data[i] += rhs.data[i] ;
	}
	return *this;
}

node& node::operator-=(const node &rhs){
	for(int i=0;i<n;i++){
		data[i] -= rhs.data[i] ;
	}
	return *this;
}

const node node::operator+(const node &rhs) const{
	return node(*this) += rhs ;
}

const node node::operator+(const int &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] += rhs ;
	}
	return result;
}

const node node::operator+(const double &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] += rhs ;
	}
	return result;
}

const node node::operator-(const int &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] += rhs ;
	}
	return result;
}

const node node::operator-(const double &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] += rhs ;
	}
	return result;
}

const node node::operator*(const double &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] *= rhs ;
	}
	return result;
}

const node node::operator*(const int &rhs) const{
	node result = *this ;
	for(int i=0;i<n;i++){
		result.data[i] *= rhs ;
	}
	return result;
}

//result of multiplying two vector is a number.
const double node::operator*(const node &rhs) const{
	double result = 0 ;
	for(int i=0;i<n;i++){
		result += (this->data[i] * rhs.data[i]) ;
	}
	return result;
}

const node node::operator-(const node &rhs) const{
	return node(*this) -= rhs ;
}

string node::toString(){
	stringstream nstr ;
	for(int i=0;i<n;i++){
		nstr << data[i] << " " ;
	}
	return nstr.str();
}
void node::dump(){
	for(int i=0;i<n;i++){
		cout << data[i] << " " ;
	}
	cout << endl ;
}

void node::addChildren(node* k){
	children.push_back(k);
}

void node::addCast(node* k){
	cast.push_back(k);
}

vector<node*> node::getChildren(){
	return children;
}

vector<node*> node::getCast(){
	return cast;
}

node* node::getChild(int i){
	return children[i];
}

node* node::getCast(int i){
	return cast[i];
}

int node::getSize(){
	return children.size();
}

void node::setParent(node* p){
	parent = p ;
}

node* node::getParent(){
	return parent;
}

bool node::isRoot(){
	return root ;
}

void node::setRoot(){
	root=true;
}

//recursively returns the nearestNode between this node and all of it's children w.r.t. q_sample
//bugfix, we cannot return or even use the q_sample as the nearest node
//todo: this is an O(n) search. Use KD-tree for faster search
//todo: utilize the time factor in the search for nearest neighbor
pair<node*, double> node::getNearestNode(node* q_sample, double* min, double* max, bool timed){
	node* nearestNode = this;
	double distance = ( timed ? q_sample->timed_distance(this, min, max) : q_sample->distance(this, min, max) );
	for(int i=0;i<children.size(); i++){
		pair<node*, double> p = children[i]->getNearestNode(q_sample, min, max, timed);
		if(p.second < distance ){
			distance = p.second ;
			nearestNode = p.first ;
		}
	}
	return make_pair(nearestNode, distance);
}

void node::save(ofstream& of){
	cout << "Printing node " << getID() << endl ;
	if(root){
		of << getID() << " -1 " <<  toString() << endl ;
	}else{
		of << getID() << " " << getParent()->getID() << " " <<  toString() << endl ;
	}
	
	for(int i=0;i<children.size();i++){
		cout << "Child " << i << endl ;
		getChild(i)->save(of);
	}
}

int node::getIndex(){
	return index;
}

void node::setIndex(int i){
	index=i;
}

void node::setCounter(int x){
	counter = x;
}

int node::getCounter(){
	return counter;
}

void node::setInputVector(vector<double> x){
	for (int i = 0; i < x.size(); i++){
		input.push_back(x[i]);
	}
}

vector<double> node::getInputVector(){
	return input;
}

double node::getInput(int i){
	return input[i];
}

void node::setJitter(int x){
	jitter = x;
}

int node::getJitter(){
	return jitter;
}